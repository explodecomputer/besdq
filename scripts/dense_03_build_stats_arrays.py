#!/usr/bin/env python3
"""Issue 03: Build beta and zscore Zarr arrays from all VCFs.

Loads the union variant table built by dense_02, then fills a (N_variants,
N_traits) float16 Zarr array for each of beta and zscore. Variants absent
from a trait are left as NaN. Applies ALID sign-flip when ALT is not the
alphabetically-first allele.

Usage:
    python dense_03_build_stats_arrays.py manifest.tsv output.besdz [compressor] [--workers N]
    compressor: none | zstd (default: none)
"""

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc

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
    manifest = read_manifest(manifest_path)
    out = Path(out_dir)
    if not out.exists():
        sys.exit(f"{out_dir} does not exist — run dense_02 first")

    print("Loading variant key index...")
    keys = np.load(str(out / "variant_keys.npy"), allow_pickle=True)
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
    p = argparse.ArgumentParser()
    p.add_argument("manifest")
    p.add_argument("out_dir")
    p.add_argument("compressor", nargs="?", default="none")
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()
    if args.compressor not in ("none", "zstd"):
        sys.exit(f"Unknown compressor '{args.compressor}': choose none or zstd")
    main(args.manifest, args.out_dir, args.compressor, args.workers)
