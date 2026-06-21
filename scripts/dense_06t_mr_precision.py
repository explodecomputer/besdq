#!/usr/bin/env python3
"""Issue 14: MR precision validation — TileDB vs original VCF.

Mirrors dense_06_mr_precision.py but reads from .besdt TileDB stores.
For each exposure/outcome pair:
  1. Extract chr1 instruments (|Z| > 5.45) from TileDB raw store.
  2. Look up beta and SE in the outcome trait from TileDB.
  3. Run IVW MR on TileDB-derived stats.
  4. Read same variants from original VCF files.
  5. Run IVW MR on VCF-derived stats.
  6. Compare: rel error < 1%, Δlog10p < 0.5.

Both raw and zstd TileDB stores are tested.

Usage:
    python dense_06t_mr_precision.py
"""

import json
from pathlib import Path

import cyvcf2
import numpy as np
import tiledb

STORE_RAW = "data/ukb-chr1.besdt"
STORE_ZSTD = "data/ukb-chr1-zstd.besdt"
STORE_RAW_ZARR = "data/ukb-chr1.besdz"  # for variant_keys and manifest path
Z_THRESH = 5.45


def load_keys(store_dir: str) -> dict:
    keys = np.load(f"{store_dir}/variant_keys.npy", allow_pickle=True)
    return {(str(r[0]), int(r[1]), str(r[2]), str(r[3])): i for i, r in enumerate(keys)}


def tiledb_col(store_dir: str, array_name: str, col: int) -> np.ndarray:
    with tiledb.open(f"{store_dir}/{array_name}", mode="r") as A:
        return A[:, col]["data"].view(np.float16).astype(np.float32)


def tiledb_rows_col(store_dir: str, array_name: str, rows: np.ndarray, col: int) -> np.ndarray:
    """Read specific rows of one column from a TileDB 2D array."""
    if len(rows) == 0:
        return np.array([], dtype=np.float32)
    with tiledb.open(f"{store_dir}/{array_name}", mode="r") as A:
        chunks = []
        i = 0
        while i < len(rows):
            run_start = rows[i]
            run_end = rows[i]
            while i + 1 < len(rows) and rows[i + 1] == rows[i] + 1:
                i += 1
                run_end = rows[i]
            # uint16 raw bits reinterpreted as float16
            chunk = A[int(run_start):int(run_end) + 1, col]["data"].view(np.float16)
            chunks.append(chunk.astype(np.float32))
            i += 1
    return np.concatenate(chunks)


def ivw_mr(beta_x, se_x, beta_y, se_y):
    w = 1.0 / se_y**2
    estimate = np.sum(w * beta_y * beta_x / se_x) / np.sum(w * beta_x**2 / se_x**2)
    se_est = np.sqrt(1.0 / np.sum(w * beta_x**2 / se_x**2))
    z = estimate / se_est
    p = 2 * (1 - _norm_cdf(abs(z)))
    return estimate, se_est, z, p


def _norm_cdf(x):
    from math import erfc, sqrt
    return 0.5 * erfc(-x / sqrt(2))


def extract_vcf_stats(vcf_path: str, positions: set) -> dict:
    result = {}
    vcf = cyvcf2.VCF(vcf_path)
    for v in vcf:
        key = (v.CHROM, v.POS)
        if key not in positions:
            continue
        es_vals = v.format("ES")
        se_vals = v.format("SE")
        if es_vals is None or se_vals is None:
            continue
        beta = float(es_vals[0][0])
        se = float(se_vals[0][0])
        if se == 0:
            continue
        alt = v.ALT[0]
        ref = v.REF
        if not (alt < ref):
            beta = -beta
        result[(v.CHROM, v.POS, min(ref, alt), max(ref, alt))] = (beta, se)
    vcf.close()
    return result


def run_pair(store_dir: str, keys_array: np.ndarray, meta: dict,
             exp_col: int, out_col: int, store_name: str, manifest: list):
    exp_id = meta["trait_ids"][exp_col]
    out_id = meta["trait_ids"][out_col]

    z_exp = tiledb_col(store_dir, "zscore", exp_col)
    inst_rows = np.where(np.abs(z_exp) > Z_THRESH)[0]
    if len(inst_rows) < 3:
        return None

    beta_x = tiledb_rows_col(store_dir, "beta", inst_rows, exp_col)
    z_x = z_exp[inst_rows].astype(np.float32)
    se_x = np.abs(beta_x / z_x)

    beta_y = tiledb_rows_col(store_dir, "beta", inst_rows, out_col)
    z_y = tiledb_rows_col(store_dir, "zscore", inst_rows, out_col)
    se_y = np.abs(beta_y / z_y)

    valid = (np.isfinite(beta_x) & np.isfinite(se_x) &
             np.isfinite(beta_y) & np.isfinite(se_y) &
             (se_x > 0) & (se_y > 0))
    if valid.sum() < 3:
        return None

    tdb_est, _, _, tdb_p = ivw_mr(
        beta_x[valid], se_x[valid], beta_y[valid], se_y[valid]
    )

    inst_keys = keys_array[inst_rows[valid]]
    positions = {(str(k[0]), int(k[1])) for k in inst_keys}

    out_vcf_path = manifest[out_col]["file_path"]
    exp_vcf_path = manifest[exp_col]["file_path"]
    vcf_out = extract_vcf_stats(out_vcf_path, positions)
    vcf_exp = extract_vcf_stats(exp_vcf_path, positions)

    beta_y_vcf, se_y_vcf, beta_x_vcf, se_x_vcf = [], [], [], []
    for k in inst_keys:
        full_key = (str(k[0]), int(k[1]), str(k[2]), str(k[3]))
        oy = vcf_out.get(full_key)
        ex = vcf_exp.get(full_key)
        if oy is None or ex is None:
            continue
        beta_y_vcf.append(oy[0])
        se_y_vcf.append(oy[1])
        beta_x_vcf.append(ex[0])
        se_x_vcf.append(ex[1])

    if len(beta_y_vcf) < 3:
        return None

    vcf_est, _, _, vcf_p = ivw_mr(
        np.array(beta_x_vcf), np.array(se_x_vcf),
        np.array(beta_y_vcf), np.array(se_y_vcf),
    )

    rel_err = abs(tdb_est - vcf_est) / (abs(vcf_est) + 1e-10)
    logp_diff = abs(np.log10(max(tdb_p, 1e-300)) - np.log10(max(vcf_p, 1e-300)))

    return {
        "store": store_name,
        "exposure": exp_id,
        "outcome": out_id,
        "n_instruments": int(valid.sum()),
        "vcf_est": float(vcf_est),
        "tdb_est": float(tdb_est),
        "rel_err": float(rel_err),
        "vcf_logp": float(np.log10(max(vcf_p, 1e-300))),
        "tdb_logp": float(np.log10(max(tdb_p, 1e-300))),
        "logp_diff": float(logp_diff),
        "pass": rel_err < 0.01 and logp_diff < 0.5,
    }


def main():
    meta = json.load(open(f"{STORE_RAW}/metadata.json"))
    manifest_path = meta["source"]
    manifest = []
    with open(manifest_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            manifest.append(dict(zip(header, line.strip().split("\t"))))

    keys_array = np.load(f"{STORE_RAW}/variant_keys.npy", allow_pickle=True)

    # Pick exposures with enough instruments using raw store
    exp_cols = []
    for col in range(100):
        z_col = tiledb_col(STORE_RAW, "zscore", col)
        if np.sum(np.abs(z_col) > Z_THRESH) >= 10:
            exp_cols.append(col)

    pairs = []
    for exp_col in exp_cols[:5]:
        out_options = [c for c in range(100) if c != exp_col][:2]
        for out_col in out_options:
            pairs.append((exp_col, out_col))

    print(f"{'Store':<16} {'Exposure':<16} {'Outcome':<16} {'N_inst':>6} "
          f"{'VCF_est':>10} {'TDB_est':>10} {'rel_err%':>9} {'logp_diff':>10} {'PASS'}")
    print("-" * 107)

    all_pass = True
    for exp_col, out_col in pairs:
        for store_dir, store_name in [
            (STORE_RAW, "raw_float16"),
            (STORE_ZSTD, "zstd_bitshuffle"),
        ]:
            r = run_pair(store_dir, keys_array, meta, exp_col, out_col, store_name, manifest)
            if r is None:
                continue
            flag = "PASS" if r["pass"] else "FAIL"
            if not r["pass"]:
                all_pass = False
            print(f"{r['store']:<16} {r['exposure']:<16} {r['outcome']:<16} "
                  f"{r['n_instruments']:>6} {r['vcf_est']:>10.4f} {r['tdb_est']:>10.4f} "
                  f"{r['rel_err']*100:>8.3f}% {r['logp_diff']:>10.3f}  {flag}")

    print()
    print("All pairs PASS" if all_pass else "Some pairs FAIL — review above")


if __name__ == "__main__":
    main()
