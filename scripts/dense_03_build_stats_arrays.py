#!/usr/bin/env python3
"""Issue 03: Build beta and zscore Zarr arrays from all VCFs.

Loads the union variant table built by dense_02, then fills a (N_variants,
N_traits) float16 Zarr array for each of beta and zscore. Variants absent
from a trait are left as NaN. Applies ALID sign-flip when ALT is not the
alphabetically-first allele.

Usage:
    python dense_03_build_stats_arrays.py manifest.tsv output.besdz [compressor]
    compressor: none | zstd (default: none)
"""

import json
import sys
from datetime import date
from pathlib import Path

import cyvcf2
import numpy as np
import zarr
from numcodecs import Blosc


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


def main(manifest_path: str, out_dir: str, compressor_name: str = "none") -> None:
    manifest = read_manifest(manifest_path)
    out = Path(out_dir)
    if not out.exists():
        sys.exit(f"{out_dir} does not exist — run dense_02 first")

    # Load variant key index and build lookup dict
    print("Loading variant key index...")
    keys = np.load(str(out / "variant_keys.npy"), allow_pickle=True)
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

    print(f"\nWriting Zarr arrays ({compressor_name} compression)...")

    if compressor_name == "zstd":
        compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    else:
        compressor = None

    chunk_shape = (1000, n_traits)
    store = zarr.open(str(out), mode="a")
    store.create_dataset(
        "beta", data=beta_arr, chunks=chunk_shape,
        compressor=compressor, dtype="float16", overwrite=True,
    )
    store.create_dataset(
        "zscore", data=zscore_arr, chunks=chunk_shape,
        compressor=compressor, dtype="float16", overwrite=True,
    )

    # Metadata JSON
    meta = {
        "trait_ids": [r["trait_id"] for r in manifest],
        "trait_names": [r["trait_name"] for r in manifest],
        "n_per_trait": [int(r["n"]) for r in manifest],
        "n_variants": n_variants,
        "n_traits": n_traits,
        "compressor": compressor_name,
        "chunk_shape": list(chunk_shape),
        "build_date": str(date.today()),
        "source": str(Path(manifest_path).resolve()),
    }
    with open(out / "metadata.json", "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"Done. Store: {out}")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(f"Usage: {sys.argv[0]} <manifest.tsv> <output.besdz> [none|zstd]")
    compressor = sys.argv[3] if len(sys.argv) == 4 else "none"
    main(sys.argv[1], sys.argv[2], compressor)
