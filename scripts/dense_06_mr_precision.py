#!/usr/bin/env python3
"""Issue 06: MR precision validation — Zarr vs original VCF.

For each exposure/outcome pair:
  1. Extract chr1 instruments (|Z| > 5.45) for the exposure from the Zarr store.
  2. Look up those variants in the outcome trait (Zarr and VCF).
  3. Run IVW MR on both; compare point estimate and log10-p.

Tolerance: <1% relative error on IVW estimate, <0.5 log10-p difference.

Usage:
    python dense_06_mr_precision.py
"""

import json
import sys
from pathlib import Path

import cyvcf2
import numpy as np
import zarr


STORE_RAW = "data/ukb-chr1.besdz"
STORE_ZSTD = "data/ukb-chr1-zstd.besdz"
Z_THRESH = 5.45  # ~ p < 5e-8


def load_keys(store_dir: str) -> dict:
    keys = np.load(f"{store_dir}/variant_keys.npy", allow_pickle=True)
    return {(str(r[0]), int(r[1]), str(r[2]), str(r[3])): i for i, r in enumerate(keys)}


def ivw_mr(beta_x, se_x, beta_y, se_y):
    """IVW MR. Returns (estimate, se, z, p)."""
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
    """Extract beta, se for variants at given (chrom, pos) positions from VCF."""
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
        # Apply ALID flip
        alt = v.ALT[0]
        ref = v.REF
        if not (alt < ref):  # flip: REF is A1
            beta = -beta
        result[(v.CHROM, v.POS, min(ref, alt), max(ref, alt))] = (beta, se)
    vcf.close()
    return result


def run_pair(store, key_to_row, meta, exp_col, out_col, store_name, manifest):
    exp_id = meta["trait_ids"][exp_col]
    out_id = meta["trait_ids"][out_col]

    # Instruments: significant variants in exposure
    z_exp = store["zscore"][:, exp_col]
    inst_rows = np.where(np.abs(z_exp) > Z_THRESH)[0]
    if len(inst_rows) < 3:
        return None

    beta_x = store["beta"][inst_rows, exp_col].astype(np.float32)
    z_x = z_exp[inst_rows].astype(np.float32)
    se_x = np.abs(beta_x / z_x)

    beta_y_zarr = store["beta"][inst_rows, out_col].astype(np.float32)
    z_y_zarr = store["zscore"][inst_rows, out_col].astype(np.float32)
    se_y_zarr = np.abs(beta_y_zarr / z_y_zarr)

    # Drop NaN in outcome
    valid = np.isfinite(beta_x) & np.isfinite(se_x) & np.isfinite(beta_y_zarr) & np.isfinite(se_y_zarr) & (se_x > 0) & (se_y_zarr > 0)
    if valid.sum() < 3:
        return None

    zarr_est, _, _, zarr_p = ivw_mr(beta_x[valid], se_x[valid], beta_y_zarr[valid], se_y_zarr[valid])

    # VCF comparison: read outcome VCF for the same variants
    keys_array = np.load(f"{STORE_RAW}/variant_keys.npy", allow_pickle=True)
    inst_keys = keys_array[inst_rows[valid]]
    positions = {(str(k[0]), int(k[1])) for k in inst_keys}

    out_vcf_path = manifest[out_col]["file_path"]
    vcf_stats = extract_vcf_stats(out_vcf_path, positions)

    # Align VCF outcome stats to instrument order
    beta_y_vcf, se_y_vcf, beta_x_vcf, se_x_vcf = [], [], [], []
    exp_vcf_path = manifest[exp_col]["file_path"]
    exp_vcf_stats = extract_vcf_stats(exp_vcf_path, positions)

    for k in inst_keys:
        full_key = (str(k[0]), int(k[1]), str(k[2]), str(k[3]))
        pos_key = (str(k[0]), int(k[1]))
        oy = vcf_stats.get(full_key)
        ex = exp_vcf_stats.get(full_key)
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

    rel_err = abs(zarr_est - vcf_est) / (abs(vcf_est) + 1e-10)
    logp_diff = abs(np.log10(max(zarr_p, 1e-300)) - np.log10(max(vcf_p, 1e-300)))

    return {
        "store": store_name,
        "exposure": exp_id,
        "outcome": out_id,
        "n_instruments": int(valid.sum()),
        "vcf_est": float(vcf_est),
        "zarr_est": float(zarr_est),
        "rel_err": float(rel_err),
        "vcf_logp": float(np.log10(max(vcf_p, 1e-300))),
        "zarr_logp": float(np.log10(max(zarr_p, 1e-300))),
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

    store_raw = zarr.open(STORE_RAW, mode="r")
    store_zstd = zarr.open(STORE_ZSTD, mode="r")
    key_to_row = load_keys(STORE_RAW)

    # Pick exposure traits with enough hits
    z_col_cache = {col: store_raw["zscore"][:, col] for col in range(100)}
    exp_cols = [col for col in range(100)
                if np.sum(np.abs(z_col_cache[col]) > Z_THRESH) >= 10]

    # Use first 5 exposure traits, each paired with 2 outcome traits
    pairs = []
    for exp_col in exp_cols[:5]:
        out_options = [c for c in range(100) if c != exp_col][:2]
        for out_col in out_options:
            pairs.append((exp_col, out_col))

    print(f"{'Store':<16} {'Exposure':<16} {'Outcome':<16} {'N_inst':>6} "
          f"{'VCF_est':>10} {'Zarr_est':>10} {'rel_err%':>9} {'logp_diff':>10} {'PASS'}")
    print("-" * 105)

    all_pass = True
    for exp_col, out_col in pairs:
        for store, store_name in [(store_raw, "raw_float16"), (store_zstd, "zstd_bitshuffle")]:
            r = run_pair(store, key_to_row, meta, exp_col, out_col, store_name, manifest)
            if r is None:
                continue
            flag = "PASS" if r["pass"] else "FAIL"
            if not r["pass"]:
                all_pass = False
            print(f"{r['store']:<16} {r['exposure']:<16} {r['outcome']:<16} "
                  f"{r['n_instruments']:>6} {r['vcf_est']:>10.4f} {r['zarr_est']:>10.4f} "
                  f"{r['rel_err']*100:>8.3f}% {r['logp_diff']:>10.3f}  {flag}")

    print()
    print("All pairs PASS" if all_pass else "Some pairs FAIL — review above")


if __name__ == "__main__":
    main()
