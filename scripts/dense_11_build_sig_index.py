#!/usr/bin/env python3
"""Issue 11: Build significance indexes for a .besdz store.

For each of three p-value thresholds (5e-8, 5e-6, 5e-4), builds a flat
uint32 array of significant variant row indices (concatenated across all
traits) and a matching int64 offset array, stored in-place in the existing
Zarr store.

Variants with NaN zscore are excluded from all indexes.

Usage:
    python scripts/dense_11_build_sig_index.py <store_path>
"""

import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

THRESHOLDS = [5e-8, 5e-6, 5e-4]
THRESHOLD_KEYS = ["5e8", "5e6", "5e4"]
# Two-sided z cutoffs: |Z| > z_cut ↔ p < threshold
Z_CUTOFFS = [5.4507, 4.8916, 3.7190]


def build_sig_index(store_path: str) -> None:
    store_dir = Path(store_path)
    if not store_dir.exists():
        sys.exit(f"Store not found: {store_dir}")

    store = zarr.open(str(store_dir), mode="r+")

    if "zscore" not in store:
        sys.exit("Store has no 'zscore' array")

    meta_path = store_dir / "metadata.json"
    meta = json.load(open(meta_path))
    n_traits = meta["n_traits"]

    print(f"Loading zscore ({store['zscore'].shape})...")
    Z = store["zscore"][:].astype(np.float32)   # (n_variants, n_traits)
    absZ = np.abs(Z)
    is_nan = np.isnan(Z)

    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)

    for thresh, key, zcut in zip(THRESHOLDS, THRESHOLD_KEYS, Z_CUTOFFS):
        print(f"  p < {thresh:.0e}  |Z| > {zcut:.4f} ...", end="  ")

        flat_parts = []
        offsets = np.zeros(n_traits + 1, dtype=np.int64)

        for t in range(n_traits):
            hits = np.where((absZ[:, t] > zcut) & ~is_nan[:, t])[0].astype(np.uint32)
            flat_parts.append(hits)
            offsets[t + 1] = offsets[t] + len(hits)

        flat_arr = np.concatenate(flat_parts) if flat_parts else np.array([], dtype=np.uint32)
        n_hits = len(flat_arr)
        print(f"{n_hits:,} hits ({n_hits / n_traits:.0f}/trait avg)")

        chunk_size = max(1, min(n_hits, 100_000))
        store.create_dataset(
            f"sig_{key}", data=flat_arr,
            chunks=(chunk_size,), compressor=compressor, dtype="uint32", overwrite=True,
        )
        store.create_dataset(
            f"sig_{key}_offsets", data=offsets,
            chunks=(n_traits + 1,), compressor=compressor, dtype="int64", overwrite=True,
        )

    meta["sig_index"] = {
        "thresholds": THRESHOLDS,
        "keys": THRESHOLD_KEYS,
        "build_date": str(date.today()),
    }
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)

    total_mb = sum(f.stat().st_size for f in store_dir.rglob("*") if f.is_file()) / 1e6
    print(f"Done. Store: {store_dir}  ({total_mb:.0f} MB)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <store_path>")
    build_sig_index(sys.argv[1])
