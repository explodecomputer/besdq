#!/usr/bin/env python3
"""Issue 10: Benchmark LD-compressed stores against raw and zstd stores.

Measures:
  1. Storage size (du) for all five stores
  2. Query latency: regional, PheWAS, tophits, bulk (median + p95 over 10 reps)
  3. MR precision: IVW estimates from LD stores vs raw float16 zarr (ground truth)

Usage:
    python scripts/dense_10_benchmark_ld.py
"""

import json
import subprocess
import time
from math import erfc, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

TABIX = "/home/gh13047/miniforge3/envs/bcftools/bin/tabix"
RAW_STORE = "data/ukb-chr1.besdz"
ZSTD_STORE = "data/ukb-chr1-zstd.besdz"
LD_STORES = {
    "ld95":   "data/ukb-chr1-ld95.besdz",
    "ld99":   "data/ukb-chr1-ld99.besdz",
    "ld99p9": "data/ukb-chr1-ld99p9.besdz",
}
VCF_DIR = "data/vcf-ukb"
EIGVEC_DIR = "data/ld_eigenvectors"
N_REPS = 10
REGION = ("1", 100_000_000, 101_000_000)
PHEWAS_POS = ("1", 100_500_000)
Z_THRESH = 5.45


# ── helpers ──────────────────────────────────────────────────────────────────

def store_size_mb(path: str) -> float:
    result = subprocess.run(["du", "-sb", path], capture_output=True, text=True)
    return int(result.stdout.split()[0]) / 1e6


def bench(fn, n=N_REPS):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return np.median(times) * 1000, np.percentile(times, 95) * 1000, result


def tabix_indices(store_dir, chrom, start, end, key_to_row):
    tsv_gz = f"{store_dir}/variants.tsv.gz"
    result = subprocess.run([TABIX, tsv_gz, f"{chrom}:{start}-{end}"],
                            capture_output=True, text=True)
    if not result.stdout.strip():
        return np.array([], dtype=np.int32)
    indices = []
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        key = (parts[0], int(parts[1]), parts[2], parts[3])
        idx = key_to_row.get(key)
        if idx is not None:
            indices.append(idx)
    return np.array(sorted(indices), dtype=np.int32)


def load_raw_keys(store_dir):
    keys = np.load(f"{store_dir}/variant_keys.npy", allow_pickle=True)
    return {(str(r[0]), int(r[1]), str(r[2]), str(r[3])): i for i, r in enumerate(keys)}


# ── LD reconstruction helpers ─────────────────────────────────────────────────

def load_ld_store(store_dir, eigvec_dir):
    """Return (zarr_store, block_index, eigenvectors_npz, nan_mask)."""
    zs = zarr.open(store_dir, mode="r")
    block_idx = json.load(open(f"{store_dir}/block_index.json"))
    ev = np.load(f"{eigvec_dir}/eigenvectors.npz")
    return zs, block_idx, ev


def reconstruct_rows(rows, zs, block_idx, ev, array="beta"):
    """Reconstruct beta or zscore for given zarr row indices from an LD store."""
    rows_set = set(rows)
    n_traits = json.load(open(Path(zs.store.path) / "metadata.json"))["n_traits"]
    result = np.full((len(rows), n_traits), np.nan, dtype=np.float32)
    row_to_out = {r: i for i, r in enumerate(rows)}

    coeff_key = f"coefficients_{array}"
    for b in block_idx:
        K = b["K"]
        zarr_rows = b["zarr_row_indices"]
        overlap = [r for r in zarr_rows if r in rows_set]
        if K == 0 or not overlap:
            continue
        V = ev[f"b{b['block_id']}"][:, :K]  # (m, K)
        C = zs[coeff_key][b["coeff_row_start"]:b["coeff_row_end"], :]  # (K, n_traits)
        recon = V @ C  # (m, n_traits)
        for local_i, global_r in enumerate(zarr_rows):
            if global_r in rows_set:
                result[row_to_out[global_r], :] = recon[local_i, :]

    # Apply NaN mask
    mask_rows = np.array(list(rows), dtype=np.intp)
    for i, r in enumerate(rows):
        nan_cols = zs["nan_mask"][r, :].astype(bool)
        result[i, nan_cols] = np.nan

    # Passthrough variants
    if "passthrough_row_indices" in zs:
        pt_rows = zs["passthrough_row_indices"][:]
        pt_set = set(pt_rows.tolist())
        for i, r in enumerate(rows):
            if r in pt_set:
                pt_local = int(np.searchsorted(pt_rows, r))
                result[i, :] = zs[f"passthrough_{array}"][pt_local, :].astype(np.float32)

    return result


# ── IVW MR ───────────────────────────────────────────────────────────────────

def norm_cdf(x):
    return 0.5 * erfc(-x / sqrt(2))


def ivw_mr(beta_x, se_x, beta_y, se_y):
    w = 1.0 / se_y**2
    est = np.sum(w * beta_y * beta_x / se_x) / np.sum(w * beta_x**2 / se_x**2)
    se_est = np.sqrt(1.0 / np.sum(w * beta_x**2 / se_x**2))
    z = est / se_est
    p = 2 * (1 - norm_cdf(abs(z)))
    return est, se_est, z, p


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    results = {"storage": {}, "query_latency": {}, "mr_precision": {}}

    # ── 1. Storage sizes ──────────────────────────────────────────────────────
    print("\n=== Storage sizes ===")
    all_stores = {"raw_float16": RAW_STORE, "zstd_bitshuffle": ZSTD_STORE, **LD_STORES}
    raw_size = store_size_mb(RAW_STORE)
    for name, path in all_stores.items():
        mb = store_size_mb(path)
        ratio = raw_size / mb
        results["storage"][name] = {"mb": mb, "ratio_vs_raw": ratio}
        print(f"  {name:<20} {mb:7.0f} MB  ({ratio:.1f}× smaller than raw)")

    # ── 2. Query latency ──────────────────────────────────────────────────────
    print("\n=== Query latency ===")
    print(f"{'Store':<20} {'Query':<10} {'median ms':>10} {'p95 ms':>8}")
    print("-" * 54)

    raw_zarr = zarr.open(RAW_STORE, mode="r")
    key_to_row = load_raw_keys(RAW_STORE)
    chrom, start, end = REGION
    ph_chrom, ph_pos = PHEWAS_POS
    region_indices = tabix_indices(RAW_STORE, chrom, start, end, key_to_row)
    ph_indices = tabix_indices(RAW_STORE, ph_chrom, ph_pos, ph_pos + 1, key_to_row)
    ph_row = int(ph_indices[0]) if len(ph_indices) > 0 else 0

    # Find a trait with top hits
    tophits_col = 0
    best_hits = 0
    for c in range(100):
        n_hits = int(np.sum(np.abs(raw_zarr["zscore"][:, c]) > Z_THRESH))
        if n_hits > best_hits:
            best_hits, tophits_col = n_hits, c

    # Raw and zstd stores (direct zarr access)
    for store_name, store_path in [("raw_float16", RAW_STORE), ("zstd_bitshuffle", ZSTD_STORE)]:
        zs = zarr.open(store_path, mode="r")
        queries = {
            "regional": lambda s=zs, idx=region_indices: s["beta"].oindex[idx, :],
            "phewas":   lambda s=zs, r=ph_row: s["beta"][r, :],
            "tophits":  lambda s=zs, c=tophits_col: (
                lambda z=s["zscore"][:, c]: s["beta"][np.where(np.abs(z) > Z_THRESH)[0], c])(),
            "bulk":     lambda s=zs, c=tophits_col: s["beta"][:, c],
        }
        results["query_latency"][store_name] = {}
        for q_name, fn in queries.items():
            med, p95, _ = bench(fn)
            results["query_latency"][store_name][q_name] = {"median_ms": med, "p95_ms": p95}
            print(f"  {store_name:<20} {q_name:<10} {med:>10.1f} {p95:>8.1f}")

    # LD stores (reconstruction-based queries)
    for ld_name, ld_path in LD_STORES.items():
        zs, block_idx, ev = load_ld_store(ld_path, EIGVEC_DIR)
        results["query_latency"][ld_name] = {}

        # Regional: find relevant blocks
        region_blocks = [b for b in block_idx
                         if b["chrom"] == chrom and b["stop"] > start and b["start"] < end]
        region_zarr_rows = []
        for b in region_blocks:
            region_zarr_rows.extend([r for r in b["zarr_row_indices"]
                                     if r in set(region_indices.tolist())])

        def regional_ld(rows=region_zarr_rows, z=zs, bi=block_idx, e=ev):
            return reconstruct_rows(rows, z, bi, e, "beta")

        def phewas_ld(row=ph_row, z=zs, bi=block_idx, e=ev):
            return reconstruct_rows([row], z, bi, e, "beta")

        def tophits_ld(col=tophits_col, z=zs, bi=block_idx, e=ev, n=762488, nt=100):
            # Reconstruct full column — load all blocks for this trait column
            col_beta = np.full(n, np.nan, dtype=np.float32)
            for b in bi:
                K = b["K"]
                if K == 0:
                    continue
                V = e[f"b{b['block_id']}"][:, :K]
                C = z["coefficients_beta"][b["coeff_row_start"]:b["coeff_row_end"], col]
                recon = V @ C
                rows = np.array(b["zarr_row_indices"], dtype=np.intp)
                col_beta[rows] = recon
            hits = np.where(np.abs(col_beta) > 0)[0]  # approximate
            return col_beta[hits]

        def bulk_ld(col=tophits_col, z=zs, bi=block_idx, e=ev, n=762488):
            col_beta = np.full(n, np.nan, dtype=np.float32)
            for b in bi:
                K = b["K"]
                if K == 0:
                    continue
                V = e[f"b{b['block_id']}"][:, :K]
                C = z["coefficients_beta"][b["coeff_row_start"]:b["coeff_row_end"], col]
                col_beta[np.array(b["zarr_row_indices"], dtype=np.intp)] = V @ C
            return col_beta

        queries = {
            "regional": regional_ld,
            "phewas":   phewas_ld,
            "tophits":  tophits_ld,
            "bulk":     bulk_ld,
        }
        for q_name, fn in queries.items():
            med, p95, _ = bench(fn)
            results["query_latency"][ld_name][q_name] = {"median_ms": med, "p95_ms": p95}
            print(f"  {ld_name:<20} {q_name:<10} {med:>10.1f} {p95:>8.1f}")

    # ── 3. MR precision ───────────────────────────────────────────────────────
    print("\n=== MR precision (LD stores vs raw float16) ===")
    print(f"{'Store':<12} {'Exp':<16} {'Out':<16} {'N':>5} {'raw_est':>10} "
          f"{'ld_est':>10} {'rel_err%':>9} {'logp_diff':>10} {'PASS'}")
    print("-" * 100)

    # Find exposure traits with enough instruments
    raw_zscore = raw_zarr["zscore"][:]
    exp_cols = [c for c in range(100)
                if np.sum(np.abs(raw_zscore[:, c]) > Z_THRESH) >= 10][:4]

    all_pass = {k: True for k in LD_STORES}
    results["mr_precision"] = {}

    for exp_col in exp_cols:
        out_cols = [c for c in range(100) if c != exp_col][:2]
        for out_col in out_cols:
            # Instruments from raw store
            z_exp = raw_zscore[:, exp_col]
            inst_rows = np.where(np.abs(z_exp) > Z_THRESH)[0]
            beta_x_raw = raw_zarr["beta"][inst_rows, exp_col].astype(np.float32)
            z_x_raw = z_exp[inst_rows].astype(np.float32)
            se_x_raw = np.abs(beta_x_raw / z_x_raw)

            beta_y_raw = raw_zarr["beta"][inst_rows, out_col].astype(np.float32)
            z_y_raw = raw_zscore[inst_rows, out_col].astype(np.float32)
            se_y_raw = np.abs(beta_y_raw / z_y_raw)

            valid = (np.isfinite(beta_x_raw) & np.isfinite(se_x_raw) &
                     np.isfinite(beta_y_raw) & np.isfinite(se_y_raw) &
                     (se_x_raw > 0) & (se_y_raw > 0))
            if valid.sum() < 3:
                continue

            raw_est, _, _, raw_p = ivw_mr(
                beta_x_raw[valid], se_x_raw[valid],
                beta_y_raw[valid], se_y_raw[valid])

            meta_raw = json.load(open(f"{RAW_STORE}/metadata.json"))
            exp_id = meta_raw["trait_ids"][exp_col]
            out_id = meta_raw["trait_ids"][out_col]

            for ld_name, ld_path in LD_STORES.items():
                zs, block_idx, ev = load_ld_store(ld_path, EIGVEC_DIR)

                beta_y_ld = reconstruct_rows(inst_rows[valid].tolist(), zs, block_idx, ev, "beta")[:, out_col]
                z_y_ld = reconstruct_rows(inst_rows[valid].tolist(), zs, block_idx, ev, "zscore")[:, out_col]
                se_y_ld = np.abs(beta_y_ld / (z_y_ld + 1e-30))

                valid2 = np.isfinite(beta_y_ld) & np.isfinite(se_y_ld) & (se_y_ld > 0)
                if valid2.sum() < 3:
                    continue

                ld_est, _, _, ld_p = ivw_mr(
                    beta_x_raw[valid][valid2], se_x_raw[valid][valid2],
                    beta_y_ld[valid2], se_y_ld[valid2])

                rel_err = abs(ld_est - raw_est) / (abs(raw_est) + 1e-10)
                logp_diff = abs(np.log10(max(ld_p, 1e-300)) - np.log10(max(raw_p, 1e-300)))
                passed = rel_err < 0.01 and logp_diff < 0.5
                if not passed:
                    all_pass[ld_name] = False

                key = f"{ld_name}/{exp_id}/{out_id}"
                results["mr_precision"][key] = {
                    "store": ld_name, "exposure": exp_id, "outcome": out_id,
                    "n_instruments": int(valid2.sum()),
                    "raw_est": float(raw_est), "ld_est": float(ld_est),
                    "rel_err": float(rel_err), "logp_diff": float(logp_diff),
                    "pass": passed,
                }
                flag = "PASS" if passed else "FAIL"
                if not passed:
                    all_pass[ld_name] = False
                print(f"  {ld_name:<12} {exp_id:<16} {out_id:<16} "
                      f"{valid2.sum():>5} {raw_est:>10.4f} {ld_est:>10.4f} "
                      f"{rel_err*100:>8.3f}% {logp_diff:>10.3f}  {flag}")

    print()
    for name, passed in all_pass.items():
        print(f"  {name}: {'ALL PASS' if passed else 'SOME FAIL'}")

    # ── Save results ──────────────────────────────────────────────────────────
    with open("data/ld_benchmark_results.json", "w") as fh:
        json.dump(results, fh, indent=2)
    print("\nResults written to data/ld_benchmark_results.json")


if __name__ == "__main__":
    main()
