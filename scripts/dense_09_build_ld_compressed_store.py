#!/usr/bin/env python3
"""Issue 09: Build LD-compressed .besdz store at a given variance-explained threshold.

Projects beta and zscore onto the top-K LD eigenvectors per LDetect block.
Variants not covered by 1KG LD reference are stored as raw float16 (passthrough).
NaN mask stored separately so reconstruction can restore untested positions.

Usage:
    python scripts/dense_09_build_ld_compressed_store.py \\
        <source_store> <eigenvector_dir> <output_store> <threshold>

    threshold: float in (0, 1), e.g. 0.95, 0.99, 0.999
"""

import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import zarr
from numcodecs import Blosc


def main(source_dir: str, eigvec_dir: str, out_dir: str, threshold: float) -> None:
    src = Path(source_dir)
    evd = Path(eigvec_dir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Load source metadata
    meta = json.load(open(src / "metadata.json"))
    n_variants = meta["n_variants"]
    n_traits = meta["n_traits"]
    print(f"Source: {n_variants:,} variants × {n_traits} traits")

    # Load block manifest and eigenvectors
    manifest = json.load(open(evd / "block_manifest.json"))
    print(f"Loading eigenvectors.npz ...")
    eigenvectors = np.load(str(evd / "eigenvectors.npz"))

    # Load full source zarr arrays into memory (float32 for projection arithmetic)
    print("Loading source zarr arrays ...")
    store = zarr.open(str(src), mode="r")
    beta_src = store["beta"][:].astype(np.float32)   # (n_variants, n_traits)
    zscore_src = store["zscore"][:].astype(np.float32)
    print(f"  Loaded. Memory: ~{beta_src.nbytes * 2 / 1e9:.2f} GB")

    # Build NaN mask from source (1 = originally NaN)
    nan_mask = np.isnan(beta_src).astype(np.uint8)  # (n_variants, n_traits)

    # Zero-fill NaN in-place for projection
    beta_src[nan_mask.astype(bool)] = 0.0
    zscore_src[np.isnan(zscore_src)] = 0.0

    # Determine K per block for this threshold
    # block_manifest stores eigenvalues; K = min eigenvectors to reach threshold
    block_index = []
    all_coeff_beta = []
    all_coeff_zscore = []
    coeff_row = 0

    passthrough_rows_all = []  # global zarr row indices for passthrough variants

    for block in manifest:
        bi = block["block_id"]
        zarr_rows = block["zarr_row_indices"]
        passthrough_rows = block["passthrough_row_indices"]
        passthrough_rows_all.extend(passthrough_rows)

        if len(zarr_rows) == 0 or block["k_stored"] == 0:
            # All variants are passthrough for this block
            block_index.append({
                "block_id": bi,
                "chrom": block["chrom"],
                "start": block["start"],
                "stop": block["stop"],
                "n_intersection": 0,
                "K": 0,
                "coeff_row_start": coeff_row,
                "coeff_row_end": coeff_row,
                "zarr_row_indices": [],
            })
            continue

        # Determine K for this threshold
        eigenvalues = np.array(block["eigenvalues"])
        k_available = block["k_stored"]
        eigenvalues_stored = eigenvalues[:k_available]
        cumvar = np.cumsum(eigenvalues_stored) / max(eigenvalues.sum(), 1e-10)
        k_use = int(np.searchsorted(cumvar, threshold)) + 1
        k_use = min(k_use, k_available)

        # Load eigenvectors for this block (shape: n_intersection × k_available)
        V = eigenvectors[f"b{bi}"][:, :k_use]  # (n_intersection, k_use)

        # Extract source stats for this block's intersection variants
        zarr_rows_arr = np.array(zarr_rows, dtype=np.intp)
        beta_block = beta_src[zarr_rows_arr, :]    # (n_intersection, n_traits)
        zscore_block = zscore_src[zarr_rows_arr, :] # (n_intersection, n_traits)

        # Project: C = V.T @ data_block → (k_use, n_traits)
        C_beta = V.T @ beta_block      # (k_use, n_traits)
        C_zscore = V.T @ zscore_block  # (k_use, n_traits)

        all_coeff_beta.append(C_beta.astype(np.float32))
        all_coeff_zscore.append(C_zscore.astype(np.float32))

        block_index.append({
            "block_id": bi,
            "chrom": block["chrom"],
            "start": block["start"],
            "stop": block["stop"],
            "n_intersection": len(zarr_rows),
            "K": k_use,
            "coeff_row_start": coeff_row,
            "coeff_row_end": coeff_row + k_use,
            "zarr_row_indices": zarr_rows,
        })
        coeff_row += k_use

        if bi % 20 == 0:
            print(f"  Block {bi:3d}: K={k_use}", end="\r")

    print(f"\n  Total coefficient rows (sum_K): {coeff_row:,}")

    # Write coefficient arrays as Zarr
    print("Writing Zarr arrays ...")
    chunk_k = min(512, max(1, coeff_row // 200))
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)

    if all_coeff_beta:
        coeff_beta_arr = np.vstack(all_coeff_beta)   # (sum_K, n_traits)
        coeff_zscore_arr = np.vstack(all_coeff_zscore)
    else:
        coeff_beta_arr = np.empty((0, n_traits), dtype=np.float32)
        coeff_zscore_arr = np.empty((0, n_traits), dtype=np.float32)

    out_store = zarr.open(str(out), mode="w")
    out_store.create_dataset(
        "coefficients_beta", data=coeff_beta_arr,
        chunks=(chunk_k, n_traits), compressor=compressor, dtype="float32", overwrite=True,
    )
    out_store.create_dataset(
        "coefficients_zscore", data=coeff_zscore_arr,
        chunks=(chunk_k, n_traits), compressor=compressor, dtype="float32", overwrite=True,
    )

    # NaN mask (bitshuffle compressed)
    out_store.create_dataset(
        "nan_mask", data=nan_mask,
        chunks=(1000, n_traits), compressor=compressor, dtype="uint8", overwrite=True,
    )

    # Passthrough arrays (raw float16 — copy from source)
    if passthrough_rows_all:
        pt_rows = np.array(sorted(set(passthrough_rows_all)), dtype=np.intp)
        pt_beta = store["beta"].oindex[pt_rows, :]    # float16
        pt_zscore = store["zscore"].oindex[pt_rows, :]
        out_store.create_dataset(
            "passthrough_beta", data=pt_beta,
            chunks=(1000, n_traits), compressor=None, dtype="float16", overwrite=True,
        )
        out_store.create_dataset(
            "passthrough_zscore", data=pt_zscore,
            chunks=(1000, n_traits), compressor=None, dtype="float16", overwrite=True,
        )
        # Store the passthrough row indices so we can align them later
        out_store.create_dataset(
            "passthrough_row_indices", data=pt_rows,
            chunks=(len(pt_rows),), compressor=compressor, dtype="int32", overwrite=True,
        )
    else:
        pt_rows = np.array([], dtype=np.int32)

    # Block index JSON
    with open(out / "block_index.json", "w") as fh:
        json.dump(block_index, fh)

    # Metadata JSON
    out_meta = {
        "trait_ids": meta["trait_ids"],
        "trait_names": meta["trait_names"],
        "n_per_trait": meta["n_per_trait"],
        "n_variants": n_variants,
        "n_traits": n_traits,
        "compression": "ld_eigenvector",
        "threshold": threshold,
        "n_blocks": len(block_index),
        "total_K": int(coeff_row),
        "n_passthrough_variants": int(len(pt_rows)),
        "eigenvector_dir": str(evd.resolve()),
        "source": str(src.resolve()),
        "build_date": str(date.today()),
    }
    with open(out / "metadata.json", "w") as fh:
        json.dump(out_meta, fh, indent=2)

    size_mb = sum(f.stat().st_size for f in out.rglob("*") if f.is_file()) / 1e6
    print(f"Done. Store: {out}  ({size_mb:.0f} MB)")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} <source_store> <eigvec_dir> <output_store> <threshold>")
    main(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))
