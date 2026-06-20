#!/usr/bin/env python3
"""Issue 05: Query latency benchmark across four query patterns and two stores.

Query patterns:
  1. regional  — all variants in a 1 Mb window, all traits
  2. phewas    — single variant (by position) across all traits
  3. tophits   — variants with |Z| > 5.45 (p<5e-8) for one trait
  4. bulk      — all variants for one trait

Usage:
    python dense_05_query_benchmark.py
"""

import subprocess
import time

import numpy as np
import pandas as pd
import zarr

TABIX = "/home/gh13047/miniforge3/envs/bcftools/bin/tabix"
STORES = {
    "raw_float16": "data/ukb-chr1.besdz",
    "zstd_bitshuffle": "data/ukb-chr1-zstd.besdz",
}
N_REPS = 10
REGION = ("1", 100_000_000, 101_000_000)  # 1 Mb window on chr1
PHEWAS_POS = ("1", 100_500_000)
TOPHITS_TRAIT = 0
Z_THRESHOLD = 5.45   # ~p < 5e-8


def tabix_row_indices(store_dir: str, chrom: str, start: int, end: int) -> np.ndarray:
    """Return row indices of variants in [start, end] via tabix."""
    tsv_gz = f"{store_dir}/variants.tsv.gz"
    result = subprocess.run(
        [TABIX, tsv_gz, f"{chrom}:{start}-{end}"],
        capture_output=True, text=True,
    )
    if not result.stdout.strip():
        return np.array([], dtype=np.int32)
    lines = result.stdout.strip().split("\n")
    # Build lookup to find row indices — we need the global row index
    # Tabix returns matching rows; we need their indices in the sorted variant table
    # Load the keys to find indices (cached after first call)
    return lines  # return raw lines for now; we'll post-process


def load_keys(store_dir: str) -> dict:
    keys = np.load(f"{store_dir}/variant_keys.npy", allow_pickle=True)
    return {(str(r[0]), int(r[1]), str(r[2]), str(r[3])): i for i, r in enumerate(keys)}


def tabix_indices(store_dir: str, chrom: str, start: int, end: int,
                  key_to_row: dict) -> np.ndarray:
    tsv_gz = f"{store_dir}/variants.tsv.gz"
    result = subprocess.run(
        [TABIX, tsv_gz, f"{chrom}:{start}-{end}"],
        capture_output=True, text=True,
    )
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


def query_regional(store, indices: np.ndarray) -> pd.DataFrame:
    if len(indices) == 0:
        return pd.DataFrame()
    beta = store["beta"].oindex[indices, :]
    zscore = store["zscore"].oindex[indices, :]
    return pd.DataFrame({"n_variants": [len(indices)], "n_traits": [beta.shape[1]]})


def query_phewas(store, row_idx: int) -> pd.DataFrame:
    beta = store["beta"][row_idx, :]
    zscore = store["zscore"][row_idx, :]
    return pd.DataFrame({"beta": beta, "zscore": zscore})


def query_tophits(store, trait_col: int, z_thresh: float) -> pd.DataFrame:
    z_col = store["zscore"][:, trait_col]
    hit_rows = np.where(np.abs(z_col) > z_thresh)[0]
    beta = store["beta"][hit_rows, trait_col]
    return pd.DataFrame({"row": hit_rows, "beta": beta, "z": z_col[hit_rows]})


def query_bulk(store, trait_col: int) -> pd.DataFrame:
    beta = store["beta"][:, trait_col]
    zscore = store["zscore"][:, trait_col]
    return pd.DataFrame({"beta": beta, "zscore": zscore})


def bench(fn, n=N_REPS):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        result = fn()
        times.append(time.perf_counter() - t0)
    return np.median(times) * 1000, np.percentile(times, 95) * 1000, result


def main():
    print(f"{'Store':<20} {'Query':<12} {'median ms':>10} {'p95 ms':>10} {'result'}")
    print("-" * 70)

    for store_name, store_dir in STORES.items():
        store = zarr.open(store_dir, mode="r")
        key_to_row = load_keys(store_dir)

        chrom, start, end = REGION
        indices = tabix_indices(store_dir, chrom, start, end, key_to_row)

        # PheWAS: find closest variant to PHEWAS_POS
        ph_chrom, ph_pos = PHEWAS_POS
        ph_indices = tabix_indices(store_dir, ph_chrom, ph_pos, ph_pos + 1, key_to_row)
        ph_row = int(ph_indices[0]) if len(ph_indices) > 0 else 0

        queries = {
            "regional":  lambda s=store, i=indices: query_regional(s, i),
            "phewas":    lambda s=store, r=ph_row: query_phewas(s, r),
            "tophits":   lambda s=store: query_tophits(s, TOPHITS_TRAIT, Z_THRESHOLD),
            "bulk":      lambda s=store: query_bulk(s, TOPHITS_TRAIT),
        }

        for q_name, fn in queries.items():
            med, p95, result = bench(fn)
            if isinstance(result, pd.DataFrame) and "n_variants" in result.columns:
                summary = f"{result['n_variants'].iloc[0]} variants"
            elif isinstance(result, pd.DataFrame):
                summary = f"{len(result)} rows"
            else:
                summary = str(result)
            print(f"{store_name:<20} {q_name:<12} {med:>10.1f} {p95:>10.1f}  {summary}")

        print()


if __name__ == "__main__":
    main()
