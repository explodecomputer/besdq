#!/usr/bin/env python3
"""Issue 05: Query latency benchmark across five query patterns and two stores.

Query patterns:
  1. regional       — all variants in a 1 Mb window, all traits
  2. phewas         — single variant (by position) across all traits
  3. tophits        — GW-significant variants for one trait (via sig index)
  4. bulk           — all variants for one trait
  5. random_lookup  — 100 random variants × 10 random traits

Usage:
    python dense_05_query_benchmark.py [--prefix ukb-chr1]
"""

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

TABIX = "/home/gh13047/miniforge3/envs/bcftools/bin/tabix"
N_REPS = 10
REGION = ("1", 100_000_000, 101_000_000)  # 1 Mb window on chr1
PHEWAS_POS = ("1", 100_500_000)
RNG = np.random.default_rng(42)
N_RANDOM_VARIANTS = 100
N_RANDOM_TRAITS = 10


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--prefix", default="ukb-chr1",
                   help="Store name prefix under data/ (default: ukb-chr1)")
    return p.parse_args()


def get_stores(prefix):
    return {
        "raw_float16":    f"data/{prefix}.besdz",
        "zstd_bitshuffle": f"data/{prefix}-zstd.besdz",
    }


def get_output_json(prefix):
    return f"data/{prefix}_zarr_benchmark.json"


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


def pick_tophits_trait(store) -> int:
    offsets = store["sig_5e8_offsets"][:]
    return int(np.argmax(np.diff(offsets)))


def query_tophits(store, trait_col: int) -> pd.DataFrame:
    offsets = store["sig_5e8_offsets"][:]
    start, end = int(offsets[trait_col]), int(offsets[trait_col + 1])
    if end <= start:
        return pd.DataFrame({"row": [], "beta": [], "z": []})
    hit_rows = store["sig_5e8"][start:end].astype(np.int64)
    beta = store["beta"].oindex[hit_rows, :][:, trait_col]
    z = store["zscore"].oindex[hit_rows, :][:, trait_col]
    return pd.DataFrame({"row": hit_rows, "beta": beta, "z": z})


def query_random_lookup(store, row_indices: np.ndarray, col_indices: np.ndarray) -> np.ndarray:
    return store["beta"].oindex[row_indices, :][:, col_indices]


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
    args = parse_args()
    STORES = get_stores(args.prefix)
    output_json = get_output_json(args.prefix)

    print(f"{'Store':<20} {'Query':<12} {'median ms':>10} {'p95 ms':>10} {'result'}")
    print("-" * 70)

    all_results = {}

    # Random indices shared across stores for a fair comparison
    first_store = next(iter(STORES.values()))
    _first = zarr.open(first_store, mode="r")
    n_variants = _first["beta"].shape[0]
    n_traits = _first["beta"].shape[1]
    rand_rows = np.sort(RNG.choice(n_variants, N_RANDOM_VARIANTS, replace=False))
    rand_cols = np.sort(RNG.choice(n_traits, N_RANDOM_TRAITS, replace=False))

    for store_name, store_dir in STORES.items():
        import os
        if not os.path.isdir(store_dir):
            print(f"  [SKIP] {store_dir} not found")
            continue
        store = zarr.open(store_dir, mode="r")
        key_to_row = load_keys(store_dir)

        chrom, start, end = REGION
        indices = tabix_indices(store_dir, chrom, start, end, key_to_row)

        ph_chrom, ph_pos = PHEWAS_POS
        ph_indices = tabix_indices(store_dir, ph_chrom, ph_pos, ph_pos + 1, key_to_row)
        ph_row = int(ph_indices[0]) if len(ph_indices) > 0 else 0

        tophits_col = pick_tophits_trait(store)

        queries = {
            "regional":      lambda s=store, i=indices: query_regional(s, i),
            "phewas":        lambda s=store, r=ph_row: query_phewas(s, r),
            "tophits":       lambda s=store, c=tophits_col: query_tophits(s, c),
            "bulk":          lambda s=store, c=tophits_col: query_bulk(s, c),
            "random_lookup": lambda s=store, ri=rand_rows, ci=rand_cols: query_random_lookup(s, ri, ci),
        }

        store_results = {}
        for q_name, fn in queries.items():
            med, p95, result = bench(fn)
            store_results[q_name] = {"median_ms": round(med, 3), "p95_ms": round(p95, 3)}
            if isinstance(result, pd.DataFrame) and "n_variants" in result.columns:
                summary = f"{result['n_variants'].iloc[0]} variants"
            elif isinstance(result, pd.DataFrame):
                summary = f"{len(result)} rows"
            elif isinstance(result, np.ndarray):
                summary = f"{result.shape} array"
            else:
                summary = str(result)
            print(f"{store_name:<20} {q_name:<16} {med:>10.1f} {p95:>10.1f}  {summary}")

        all_results[store_name] = store_results
        print()

    with open(output_json, "w") as fh:
        json.dump(all_results, fh, indent=2)
    print(f"Results written to {output_json}")


if __name__ == "__main__":
    main()
