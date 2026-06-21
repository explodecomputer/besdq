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

import json
import shutil
import sys
from datetime import date
from pathlib import Path

import cyvcf2
import numpy as np
import tiledb


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


def main(manifest_path: str, out_dir: str, compressor_name: str = "none") -> None:
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
    print(f"  {n_variants:,} variants × {n_traits} traits")

    # Allocate in-memory arrays (NaN = not tested)
    beta_arr = np.full((n_variants, n_traits), np.nan, dtype=np.float16)
    zscore_arr = np.full((n_variants, n_traits), np.nan, dtype=np.float16)

    # Fill arrays trait by trait
    for col_idx, row in enumerate(manifest):
        trait_id = row["trait_id"]
        vcf_path = row["file_path"]
        print(f"[{col_idx + 1}/{n_traits}] {trait_id}", end="\r", flush=True)

        vcf = cyvcf2.VCF(vcf_path)
        for variant in vcf:
            alts = variant.ALT
            if len(alts) != 1:
                continue
            ref = variant.REF
            alt = alts[0]
            chrom = variant.CHROM
            pos = variant.POS

            a1, a2, flip = alid_flip(ref, alt)
            key = (chrom, pos, a1, a2)
            row_idx = key_to_row.get(key)
            if row_idx is None:
                continue

            es_vals = variant.format("ES")
            se_vals = variant.format("SE")
            if es_vals is None or se_vals is None:
                continue

            beta = float(es_vals[0][0])
            se = float(se_vals[0][0])
            if se == 0.0:
                continue
            z = beta / se

            if flip:
                beta = -beta
                z = -z

            beta_arr[row_idx, col_idx] = beta
            zscore_arr[row_idx, col_idx] = z

        vcf.close()

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
    if len(sys.argv) not in (3, 4):
        sys.exit(f"Usage: {sys.argv[0]} <manifest.tsv> <output.besdt> [none|zstd]")
    compressor = sys.argv[3] if len(sys.argv) == 4 else "none"
    main(sys.argv[1], sys.argv[2], compressor)
