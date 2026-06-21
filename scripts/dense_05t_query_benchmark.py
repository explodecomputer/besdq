#!/usr/bin/env python3
"""Issue 13: Query latency benchmark for TileDB stores.

Mirrors dense_05_query_benchmark.py but reads from .besdt TileDB stores.
Five query patterns:
  1. regional       — all variants in a 1 Mb window, all traits
  2. phewas         — single variant across all traits
  3. tophits        — GW-significant variants for one trait (via sig index)
  4. bulk           — all variants for one trait
  5. random_lookup  — 100 random variants × 10 random traits

Non-contiguous row reads (regional query):
  TileDB has no oindex equivalent. Row indices from tabix are grouped into
  contiguous runs and one slice per run is issued, then results concatenated.
  This faithfully reflects TileDB's actual access pattern.

Output: data/tiledb_query_benchmark_results.json

Usage:
    python dense_05t_query_benchmark.py
"""

import json
import subprocess
import time

import numpy as np
import pandas as pd
import tiledb

TABIX = "/home/gh13047/miniforge3/envs/bcftools/bin/tabix"
STORES = {
    "raw_float16": "data/ukb-chr1.besdt",
    "zstd_bitshuffle": "data/ukb-chr1-zstd.besdt",
}
N_REPS = 10
REGION = ("1", 100_000_000, 101_000_000)
PHEWAS_POS = ("1", 100_500_000)
RNG = np.random.default_rng(42)
N_RANDOM_VARIANTS = 100
N_RANDOM_TRAITS = 10
OUTPUT_JSON = "data/tiledb_query_benchmark_results.json"


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


def tiledb_read_rows(uri: str, indices: np.ndarray) -> np.ndarray:
    """Read non-contiguous rows from a TileDB 2D array by grouping into runs.

    TileDB multi_index on dense arrays requires contiguous ranges. We sort
    indices, identify contiguous runs, issue one slice per run, concatenate.
    This is the honest implementation — it shows TileDB's real access pattern
    rather than hiding a scatter-gather behind an abstraction.
    """
    if len(indices) == 0:
        return np.empty((0, 0), dtype=np.float16)

    with tiledb.open(uri, mode="r") as A:
        n_cols = A.schema.domain.dim(1).domain[1] + 1
        chunks = []
        i = 0
        while i < len(indices):
            run_start = indices[i]
            run_end = indices[i]
            while i + 1 < len(indices) and indices[i + 1] == indices[i] + 1:
                i += 1
                run_end = indices[i]
            # uint16 raw bits reinterpreted as float16 (TileDB has no native float16)
            chunk = A[int(run_start):int(run_end) + 1, :]["data"].view(np.float16)
            chunks.append(chunk)
            i += 1
    return np.vstack(chunks) if chunks else np.empty((0, n_cols), dtype=np.float16)


def query_regional(store_dir: str, indices: np.ndarray) -> dict:
    if len(indices) == 0:
        return {"n_variants": 0}
    beta = tiledb_read_rows(f"{store_dir}/beta", indices)
    zscore = tiledb_read_rows(f"{store_dir}/zscore", indices)
    return {"n_variants": len(indices), "n_traits": beta.shape[1]}


def query_phewas(store_dir: str, row_idx: int) -> pd.DataFrame:
    with tiledb.open(f"{store_dir}/beta", mode="r") as A:
        beta = A[row_idx, :]["data"].view(np.float16)
    with tiledb.open(f"{store_dir}/zscore", mode="r") as A:
        zscore = A[row_idx, :]["data"].view(np.float16)
    return pd.DataFrame({"beta": beta, "zscore": zscore})


def pick_tophits_trait(store_dir: str) -> int:
    with tiledb.open(f"{store_dir}/sig_5e8_offsets", mode="r") as A:
        offsets = A[:]["data"]
    return int(np.argmax(np.diff(offsets)))


def query_tophits(store_dir: str, trait_col: int) -> pd.DataFrame:
    with tiledb.open(f"{store_dir}/sig_5e8_offsets", mode="r") as A:
        offsets = A[:]["data"]
    start, end = int(offsets[trait_col]), int(offsets[trait_col + 1])
    if end <= start:
        return pd.DataFrame({"row": [], "beta": [], "z": []})
    with tiledb.open(f"{store_dir}/sig_5e8", mode="r") as A:
        hit_rows = A[start:end]["data"].astype(np.int64)
    beta = tiledb_read_rows(f"{store_dir}/beta", hit_rows)[:, trait_col]
    z_vals = tiledb_read_rows(f"{store_dir}/zscore", hit_rows)[:, trait_col]
    return pd.DataFrame({"row": hit_rows, "beta": beta, "z": z_vals})


def query_random_lookup(store_dir: str, row_indices: np.ndarray, col_indices: np.ndarray) -> np.ndarray:
    result = tiledb_read_rows(f"{store_dir}/beta", row_indices)
    return result[:, col_indices]


def query_bulk(store_dir: str, trait_col: int) -> pd.DataFrame:
    with tiledb.open(f"{store_dir}/beta", mode="r") as A:
        beta = A[:, trait_col]["data"].view(np.float16)
    with tiledb.open(f"{store_dir}/zscore", mode="r") as A:
        zscore = A[:, trait_col]["data"].view(np.float16)
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

    all_results = {}

    # Random indices fixed across stores for a fair comparison (same seed as dense_05)
    first_store = next(iter(STORES.values()))
    meta = json.load(open(f"{first_store}/metadata.json"))
    n_variants, n_traits = meta["n_variants"], meta["n_traits"]
    rand_rows = np.sort(RNG.choice(n_variants, N_RANDOM_VARIANTS, replace=False))
    rand_cols = np.sort(RNG.choice(n_traits, N_RANDOM_TRAITS, replace=False))

    for store_name, store_dir in STORES.items():
        import os
        if not os.path.isdir(store_dir):
            print(f"  [SKIP] {store_dir} not found")
            continue

        key_to_row = load_keys(store_dir)
        chrom, start, end = REGION
        indices = tabix_indices(store_dir, chrom, start, end, key_to_row)

        ph_chrom, ph_pos = PHEWAS_POS
        ph_indices = tabix_indices(store_dir, ph_chrom, ph_pos, ph_pos + 1, key_to_row)
        ph_row = int(ph_indices[0]) if len(ph_indices) > 0 else 0

        tophits_col = pick_tophits_trait(store_dir)

        queries = {
            "regional":      lambda d=store_dir, i=indices: query_regional(d, i),
            "phewas":        lambda d=store_dir, r=ph_row: query_phewas(d, r),
            "tophits":       lambda d=store_dir, c=tophits_col: query_tophits(d, c),
            "bulk":          lambda d=store_dir, c=tophits_col: query_bulk(d, c),
            "random_lookup": lambda d=store_dir, ri=rand_rows, ci=rand_cols: query_random_lookup(d, ri, ci),
        }

        store_results = {}
        for q_name, fn in queries.items():
            med, p95, result = bench(fn)
            store_results[q_name] = {"median_ms": round(med, 3), "p95_ms": round(p95, 3)}
            if isinstance(result, dict) and "n_variants" in result:
                summary = f"{result['n_variants']} variants"
            elif isinstance(result, pd.DataFrame):
                summary = f"{len(result)} rows"
            elif isinstance(result, np.ndarray):
                summary = f"{result.shape} array"
            else:
                summary = str(result)
            print(f"{store_name:<20} {q_name:<16} {med:>10.1f} {p95:>10.1f}  {summary}")

        all_results[store_name] = store_results
        print()

    with open(OUTPUT_JSON, "w") as fh:
        json.dump(all_results, fh, indent=2)
    print(f"Results written to {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
