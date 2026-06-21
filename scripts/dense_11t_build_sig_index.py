#!/usr/bin/env python3
"""Issue 12: Build significance indexes for a .besdt TileDB store.

Mirrors dense_11_build_sig_index.py but reads the zscore TileDB array and
writes the index arrays as TileDB 1D dense arrays in-place.

For each of three p-value thresholds (5e-8, 5e-6, 5e-4):
  sig_{key}/          — TileDB dense 1D uint32 (flat hit row indices)
  sig_{key}_offsets/  — TileDB dense 1D int64  (per-trait start/end offsets)

Usage:
    python scripts/dense_11t_build_sig_index.py <store_path>
"""

import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import tiledb

THRESHOLDS = [5e-8, 5e-6, 5e-4]
THRESHOLD_KEYS = ["5e8", "5e6", "5e4"]
Z_CUTOFFS = [5.4507, 4.8916, 3.7190]


def create_1d_tiledb_array(uri: str, n: int, dtype: np.dtype) -> None:
    filters = tiledb.FilterList([
        tiledb.BitShuffleFilter(),
        tiledb.ZstdFilter(level=3),
    ])
    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim("rows", domain=(0, max(n - 1, 0)), tile=max(n, 1), dtype=np.int64),
        ),
        attrs=[tiledb.Attr(name="data", dtype=dtype, filters=filters)],
        cell_order="row-major",
        tile_order="row-major",
    )
    tiledb.Array.create(uri, schema)


def build_sig_index(store_path: str) -> None:
    store_dir = Path(store_path)
    if not store_dir.exists():
        sys.exit(f"Store not found: {store_dir}")

    zscore_uri = str(store_dir / "zscore")
    if tiledb.object_type(zscore_uri) != "array":
        sys.exit("Store has no 'zscore' TileDB array")

    meta_path = store_dir / "metadata.json"
    meta = json.load(open(meta_path))
    n_traits = meta["n_traits"]

    print(f"Loading zscore array...")
    with tiledb.open(zscore_uri, mode="r") as A:
        # uint16 raw bits reinterpreted as float16 (TileDB has no native float16)
        Z = A[:, :]["data"].view(np.float16).astype(np.float32)  # (n_variants, n_traits)

    absZ = np.abs(Z)
    is_nan = np.isnan(Z)

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

        for array_name, data in [
            (f"sig_{key}", flat_arr),
            (f"sig_{key}_offsets", offsets),
        ]:
            uri = str(store_dir / array_name)
            if tiledb.object_type(uri) == "array":
                tiledb.remove(uri)
            n = len(data)
            create_1d_tiledb_array(uri, n, data.dtype)
            if n > 0:
                with tiledb.open(uri, mode="w") as A:
                    A[:] = {"data": data}

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
