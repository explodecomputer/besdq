#!/usr/bin/env python3
"""Chunk configuration benchmark for .besdz stores.

Rechunks the existing 1000-trait besdz store into each candidate
chunk shape, then measures query latency for all five query patterns.

Usage:
    python scripts/dense_chunk_benchmark.py [--source data/ukb-chr1-1000.besdz]
"""

import argparse
import json
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import zarr

TABIX = "/home/gh13047/miniforge3/envs/bcftools/bin/tabix"
N_REPS = 10
REGION = ("1", 100_000_000, 101_000_000)
PHEWAS_POS = ("1", 100_500_000)
RNG = np.random.default_rng(42)
N_RANDOM_VARIANTS = 100
N_RANDOM_TRAITS = 10

CHUNK_CONFIGS = [
    (512,   1000),
    (1000,  1000),
    (4096,  128),
    (8192,  64),
    (32768, 16),
]


# ---------------------------------------------------------------------------
# Store creation
# ---------------------------------------------------------------------------

def rechunk_store(src_dir: Path, dst_dir: Path, chunk_shape: tuple) -> None:
    """Copy a besdz store with new chunk shapes for beta and zscore."""
    src = zarr.open(str(src_dir), mode="r")
    dst = zarr.open(str(dst_dir), mode="w")

    # Rechunk main arrays
    for name in ("beta", "zscore"):
        arr = src[name][:]          # load fully into RAM once
        dst.create_dataset(
            name, data=arr,
            chunks=chunk_shape,
            compressor=None,
            dtype="float16",
            overwrite=True,
        )

    # Copy sig-index and metadata verbatim (no rechunking needed)
    for name in src.array_keys():
        if name not in ("beta", "zscore"):
            zarr.copy(src[name], dst, name=name)

    # Copy sidecar files
    for fname in ("metadata.json", "variant_keys.npy", "variants.tsv.gz", "variants.tsv.gz.tbi"):
        src_f = src_dir / fname
        if src_f.exists():
            shutil.copy2(str(src_f), str(dst_dir / fname))


# ---------------------------------------------------------------------------
# Query helpers (mirrors dense_05_query_benchmark.py)
# ---------------------------------------------------------------------------

def load_keys(store_dir: Path) -> dict:
    keys = np.load(str(store_dir / "variant_keys.npy"), allow_pickle=True)
    return {(str(r[0]), int(r[1]), str(r[2]), str(r[3])): i for i, r in enumerate(keys)}


def tabix_indices(store_dir: Path, chrom: str, start: int, end: int, key_to_row: dict) -> np.ndarray:
    tsv_gz = str(store_dir / "variants.tsv.gz")
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


def pick_tophits_trait(store) -> int:
    offsets = store["sig_5e8_offsets"][:]
    return int(np.argmax(np.diff(offsets)))


def bench(fn, n=N_REPS):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return np.median(times) * 1000, np.percentile(times, 95) * 1000


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default="data/ukb-chr1-1000.besdz",
                        help="Source besdz store (default: data/ukb-chr1-1000.besdz)")
    args = parser.parse_args()

    src_dir = Path(args.source)
    if not src_dir.exists():
        raise SystemExit(f"Source store not found: {src_dir}")

    # Pre-compute shared indices from the source store (same for all configs)
    print("Pre-computing query indices from source store...")
    key_to_row = load_keys(src_dir)
    chrom, start, end = REGION
    regional_indices = tabix_indices(src_dir, chrom, start, end, key_to_row)
    print(f"  Regional query: {len(regional_indices)} variants in {chrom}:{start:,}-{end:,}")

    ph_chrom, ph_pos = PHEWAS_POS
    ph_indices = tabix_indices(src_dir, ph_chrom, ph_pos, ph_pos + 1, key_to_row)
    ph_row = int(ph_indices[0]) if len(ph_indices) > 0 else 0
    print(f"  PheWAS row: {ph_row}")

    src_store = zarr.open(str(src_dir), mode="r")
    n_variants, n_traits = src_store["beta"].shape
    tophits_col = pick_tophits_trait(src_store)
    print(f"  Top-hits trait col: {tophits_col}")

    rand_rows = np.sort(RNG.choice(n_variants, N_RANDOM_VARIANTS, replace=False))
    rand_cols = np.sort(RNG.choice(n_traits, N_RANDOM_TRAITS, replace=False))

    # -----------------------------------------------------------------------
    # Benchmark each chunk config
    # -----------------------------------------------------------------------
    QUERY_NAMES = ["regional", "phewas", "tophits", "bulk", "random_lookup"]
    results = {}

    col_w = 18
    hdr = f"{'chunk':>{col_w}}" + "".join(f"  {q:>14}" for q in QUERY_NAMES)
    sep = "-" * len(hdr)
    sub = " " * col_w + "".join(f"  {'med/p95 ms':>14}" for _ in QUERY_NAMES)
    print(f"\n{hdr}")
    print(sub)
    print(sep)

    with tempfile.TemporaryDirectory() as tmproot:
        for chunk_shape in CHUNK_CONFIGS:
            label = f"({chunk_shape[0]}, {chunk_shape[1]})"
            is_current = (src_store["beta"].chunks == chunk_shape)

            dst_dir = Path(tmproot) / label.replace(" ", "").replace("(", "").replace(")", "").replace(",", "x")
            rechunk_store(src_dir, dst_dir, chunk_shape)
            store = zarr.open(str(dst_dir), mode="r")

            queries = {
                "regional":      lambda s=store, i=regional_indices: s["beta"].oindex[i, :],
                "phewas":        lambda s=store, r=ph_row:            s["beta"][r, :],
                "tophits":       lambda s=store, c=tophits_col: (
                                     lambda off=s["sig_5e8_offsets"][:]: s["sig_5e8"][int(off[c]):int(off[c+1])]
                                 )(),
                "bulk":          lambda s=store, c=tophits_col:        s["beta"][:, c],
                "random_lookup": lambda s=store, ri=rand_rows, ci=rand_cols: s["beta"].oindex[ri, :][:, ci],
            }

            row_results = {}
            cells = []
            for q_name in QUERY_NAMES:
                med, p95 = bench(queries[q_name])
                row_results[q_name] = {"median_ms": round(med, 3), "p95_ms": round(p95, 3)}
                cells.append(f"{med:6.1f}/{p95:6.1f}")

            tag = " ← current" if is_current else ""
            print(f"{label + tag:>{col_w}}" + "".join(f"  {c:>14}" for c in cells))
            results[label] = row_results

    # -----------------------------------------------------------------------
    # Save JSON
    # -----------------------------------------------------------------------
    out_path = src_dir.parent / "chunk_benchmark_results.json"
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nResults saved → {out_path}")


if __name__ == "__main__":
    main()
