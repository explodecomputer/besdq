#!/usr/bin/env python3
"""Issue 11: Build beta and zscore TileDB arrays from all VCFs.

Mirrors dense_03_build_stats_arrays.py but writes a .besdt TileDB store
instead of a .besdz Zarr store. Array layout:

  <store>/beta/    — TileDB dense 2D float16 (N_variants × N_traits)
  <store>/zscore/  — TileDB dense 2D float16 (N_variants × N_traits)
  <store>/metadata.json
  <store>/variant_keys.npy   (copied from .besdz sibling)
  <store>/variants.tsv.gz    (copied)
  <store>/variants.tsv.gz.tbi (copied)

Usage:
    python dense_03t_build_stats_arrays.py manifest.tsv output.besdt [compressor]
    compressor: none | zstd (default: none)
"""

import argparse
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

import numpy as np
import tiledb

BCFTOOLS = "/home/gh13047/miniforge3/envs/bcftools/bin/bcftools"


def alid_flip(ref: str, alt: str) -> tuple[str, str, bool]:
    """Return (A1, A2, flip). flip=True means negate beta/z."""
    if alt < ref:
        return alt, ref, False
    return ref, alt, True


def read_manifest(manifest_path: str) -> list[dict]:
    rows = []
    with open(manifest_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def create_tiledb_array(
    uri: str,
    n_rows: int,
    n_cols: int,
    compressor_name: str,
) -> None:
    """Create an empty TileDB dense 2D array storing float16 bits as uint16.

    TileDB has no native float16 dtype. We store raw float16 bit patterns as
    uint16 (identical byte layout) and reinterpret to float16 on read.
    """
    if compressor_name == "zstd":
        filters = tiledb.FilterList([
            tiledb.BitShuffleFilter(),
            tiledb.ZstdFilter(level=3),
        ])
    else:
        filters = tiledb.FilterList([])

    schema = tiledb.ArraySchema(
        domain=tiledb.Domain(
            tiledb.Dim("rows", domain=(0, n_rows - 1), tile=1000, dtype=np.int64),
            tiledb.Dim("cols", domain=(0, n_cols - 1), tile=n_cols, dtype=np.int64),
        ),
        attrs=[tiledb.Attr(name="data", dtype=np.uint16, filters=filters)],
        cell_order="row-major",
        tile_order="row-major",
    )
    tiledb.Array.create(uri, schema)


def write_tiledb_array(uri: str, data: np.ndarray) -> None:
    """Write a 2D float16 numpy array into a TileDB array (stored as uint16 bits)."""
    with tiledb.open(uri, mode="w") as A:
        A[:, :] = {"data": data.view(np.uint16)}


def _fill_col(col_idx: int, vcf_path: str, key_to_row: dict,
              beta_arr: np.ndarray, zscore_arr: np.ndarray) -> int:
    """Fill one column of beta_arr/zscore_arr from a single VCF."""
    cmd = [BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%ES\t%SE]\n", vcf_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        chrom, pos_str, ref, alt = parts[0], parts[1], parts[2], parts[3]
        if "," in alt:
            continue
        try:
            beta = float(parts[4])
            se = float(parts[5])
        except ValueError:
            continue
        if se == 0.0:
            continue

        a1, a2, flip = alid_flip(ref, alt)
        key = (chrom, int(pos_str), a1, a2)
        row_idx = key_to_row.get(key)
        if row_idx is None:
            continue

        z = beta / se
        if flip:
            beta = -beta
            z = -z

        beta_arr[row_idx, col_idx] = beta
        zscore_arr[row_idx, col_idx] = z

    proc.wait()
    if proc.returncode != 0:
        print(f"\nWarning: bcftools failed (exit {proc.returncode}) for {vcf_path}", file=sys.stderr)
    return col_idx


def main(manifest_path: str, out_dir: str, compressor_name: str = "none", workers: int = 1) -> None:
    if compressor_name not in ("none", "zstd"):
        sys.exit(f"Unknown compressor '{compressor_name}': choose none or zstd")

    manifest = read_manifest(manifest_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Locate .besdz sibling to copy variant lookup files from
    besdz_sibling = out.with_suffix(".besdz")
    if not besdz_sibling.exists():
        # Try stripping -zstd from name to find the raw .besdt sibling's .besdz
        raw_name = out.stem.replace("-zstd", "")
        besdz_sibling = out.parent / f"{raw_name}.besdz"
    if not besdz_sibling.exists():
        sys.exit(f"Could not find .besdz sibling for variant_keys; expected {besdz_sibling}")

    print("Loading variant key index...")
    keys = np.load(str(besdz_sibling / "variant_keys.npy"), allow_pickle=True)
    key_to_row: dict[tuple, int] = {
        (str(row[0]), int(row[1]), str(row[2]), str(row[3])): i
        for i, row in enumerate(keys)
    }
    n_variants = len(keys)
    n_traits = len(manifest)
    print(f"  {n_variants:,} variants × {n_traits} traits  (workers={workers})")

    beta_arr = np.full((n_variants, n_traits), np.nan, dtype=np.float16)
    zscore_arr = np.full((n_variants, n_traits), np.nan, dtype=np.float16)

    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_fill_col, col_idx, row["file_path"], key_to_row, beta_arr, zscore_arr): col_idx
            for col_idx, row in enumerate(manifest)
        }
        for fut in as_completed(futures):
            fut.result()
            done += 1
            print(f"  {done}/{n_traits} traits done", end="\r", flush=True)

    print(f"\nCreating TileDB arrays ({compressor_name} compression)...")

    beta_uri = str(out / "beta")
    zscore_uri = str(out / "zscore")

    for uri in (beta_uri, zscore_uri):
        if tiledb.object_type(uri) == "array":
            tiledb.remove(uri)

    create_tiledb_array(beta_uri, n_variants, n_traits, compressor_name)
    create_tiledb_array(zscore_uri, n_variants, n_traits, compressor_name)

    print("  Writing beta...")
    write_tiledb_array(beta_uri, beta_arr)
    print("  Writing zscore...")
    write_tiledb_array(zscore_uri, zscore_arr)

    # Copy variant lookup files from .besdz sibling
    for fname in ("variant_keys.npy", "variants.tsv.gz", "variants.tsv.gz.tbi"):
        src = besdz_sibling / fname
        dst = out / fname
        if src.exists() and not dst.exists():
            shutil.copy2(str(src), str(dst))

    # Metadata JSON
    chunk_shape = [1000, n_traits]
    meta = {
        "trait_ids": [r["trait_id"] for r in manifest],
        "trait_names": [r["trait_name"] for r in manifest],
        "n_per_trait": [int(r["n"]) for r in manifest],
        "n_variants": n_variants,
        "n_traits": n_traits,
        "compressor": compressor_name,
        "chunk_shape": chunk_shape,
        "build_date": str(date.today()),
        "source": str(Path(manifest_path).resolve()),
    }
    with open(out / "metadata.json", "w") as fh:
        json.dump(meta, fh, indent=2)

    total_mb = sum(f.stat().st_size for f in out.rglob("*") if f.is_file()) / 1e6
    print(f"Done. Store: {out}  ({total_mb:.0f} MB)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("out_dir")
    p.add_argument("compressor", nargs="?", default="none")
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()
    if args.compressor not in ("none", "zstd"):
        sys.exit(f"Unknown compressor '{args.compressor}': choose none or zstd")
    main(args.manifest, args.out_dir, args.compressor, args.workers)
