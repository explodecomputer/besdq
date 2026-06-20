#!/usr/bin/env python3
"""Issue 08: Compute LD eigenvectors for LDetect chr1 blocks from 1KG EUR.

For each block: load 1KG EUR genotypes (bed-reader), filter MAF >= 0.01, intersect
with UKB zarr variants, compute Pearson R, eigendecompose. Stores V truncated to
99.9% variance explained (or rank of genotype matrix, whichever is smaller).

Usage:
    python scripts/dense_08_build_ld_eigenvectors.py
"""

import bisect
import json
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

PLINK_PREFIX = "/local-scratch/data/1kg/EUR"
LDETECT = "data/ldetect_chr1.bed"
ZARR_STORE = "data/ukb-chr1.besdz"
OUT_DIR = Path("data/ld_eigenvectors")
MAF_MIN = 0.01
MIN_VARIANTS = 5       # skip blocks with fewer intersection variants
VAR_STORE_THRESH = 0.999  # store eigenvectors up to this variance explained


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    blocks = pd.read_csv(LDETECT, sep="\t")
    print(f"LDetect chr1 blocks: {len(blocks)}")

    # Load UKB zarr variant keys and build chr1 lookup
    keys = np.load(f"{ZARR_STORE}/variant_keys.npy", allow_pickle=True)
    chr1_zarr = [
        (int(r[1]), str(r[2]), str(r[3]), i)
        for i, r in enumerate(keys) if str(r[0]) == "1"
    ]
    chr1_zarr.sort()
    chr1_zarr_pos = [x[0] for x in chr1_zarr]
    # (pos, a1, a2) -> zarr_row_index
    zarr_pos_a1_a2 = {(x[0], x[1], x[2]): x[3] for x in chr1_zarr}
    print(f"UKB zarr chr1 variants: {len(chr1_zarr):,}")

    # Load 1KG EUR bim — chr1 only, keeping original row indices for bed-reader
    print("Loading 1KG EUR bim (chr1)...")
    bim_full = pd.read_csv(
        f"{PLINK_PREFIX}.bim", sep="\t",
        names=["chr", "rsid", "cm", "pos", "a1", "a2"],
        dtype={"chr": str, "a1": str, "a2": str, "pos": np.int32},
    )
    bim_full["orig_idx"] = np.arange(len(bim_full), dtype=np.int32)
    bim_chr1 = bim_full[bim_full["chr"] == "1"].copy().reset_index(drop=True)
    print(f"  1KG EUR chr1 variants: {len(bim_chr1):,}")

    # Compute ALID keys and flip flags (vectorized)
    a1 = bim_chr1["a1"].values
    a2 = bim_chr1["a2"].values
    # Skip multi-char alleles (indels)
    snp_mask = np.array([len(x) == 1 and len(y) == 1 for x, y in zip(a1, a2)])
    bim_chr1 = bim_chr1[snp_mask].copy().reset_index(drop=True)
    a1 = bim_chr1["a1"].values
    a2 = bim_chr1["a2"].values

    flip = a1 > a2  # True -> a2 is ALID A1, need to flip dosage
    bim_chr1["alid_a1"] = np.where(~flip, a1, a2)
    bim_chr1["alid_a2"] = np.where(~flip, a2, a1)
    bim_chr1["flip"] = flip

    # Pre-sort by position for block lookups
    bim_chr1 = bim_chr1.sort_values("pos").reset_index(drop=True)
    bim_pos = bim_chr1["pos"].values

    bed = open_bed(f"{PLINK_PREFIX}.bed")
    n_samples = bed.iid_count
    print(f"  1KG EUR samples: {n_samples}")

    block_manifest = []
    eigenvectors = {}
    total_compressed = 0
    total_passthrough = 0

    for bi, block in blocks.iterrows():
        chrom = str(block["chr"])
        start = int(block["start"])
        stop = int(block["stop"])

        # UKB zarr variants in this block
        lo_z = bisect.bisect_right(chr1_zarr_pos, start)
        hi_z = bisect.bisect_right(chr1_zarr_pos, stop)
        block_zarr = chr1_zarr[lo_z:hi_z]

        # 1KG variants in this block
        lo_b = bisect.bisect_right(bim_pos, start)
        hi_b = bisect.bisect_right(bim_pos, stop)
        block_bim = bim_chr1.iloc[lo_b:hi_b]

        # Find intersection
        zarr_rows = []
        snp_orig_indices = []
        flip_flags = []
        for _, brow in block_bim.iterrows():
            key = (int(brow["pos"]), str(brow["alid_a1"]), str(brow["alid_a2"]))
            zrow = zarr_pos_a1_a2.get(key)
            if zrow is not None:
                zarr_rows.append(zrow)
                snp_orig_indices.append(int(brow["orig_idx"]))
                flip_flags.append(bool(brow["flip"]))

        all_zarr_rows = [x[3] for x in block_zarr]
        zarr_rows_set = set(zarr_rows)
        passthrough_rows = [r for r in all_zarr_rows if r not in zarr_rows_set]

        if len(zarr_rows) < MIN_VARIANTS:
            # No eigenvectors for this block; all zarr variants go to passthrough
            block_manifest.append({
                "block_id": bi,
                "chrom": chrom,
                "start": start,
                "stop": stop,
                "n_1kg_variants": len(block_bim),
                "n_zarr_variants": len(all_zarr_rows),
                "n_intersection_variants": 0,
                "n_passthrough_variants": len(all_zarr_rows),
                "zarr_row_indices": [],
                "passthrough_row_indices": all_zarr_rows,
                "eigenvalues": [],
                "k_stored": 0,
            })
            total_passthrough += len(all_zarr_rows)
            continue

        # Load genotypes for intersection SNPs
        snp_idx_arr = np.array(snp_orig_indices, dtype=np.intp)
        G = bed.read(np.s_[:, snp_idx_arr]).astype(np.float32)  # (n_samples, n_snps)

        # Mean-impute missing genotypes (NaN)
        col_means = np.nanmean(G, axis=0)
        nan_mask = np.isnan(G)
        G[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

        # MAF filter
        freq = G.mean(axis=0) / 2.0
        maf = np.minimum(freq, 1.0 - freq)
        keep = maf >= MAF_MIN
        if keep.sum() < MIN_VARIANTS:
            block_manifest.append({
                "block_id": bi, "chrom": chrom, "start": start, "stop": stop,
                "n_1kg_variants": len(block_bim), "n_zarr_variants": len(all_zarr_rows),
                "n_intersection_variants": 0, "n_passthrough_variants": len(all_zarr_rows),
                "zarr_row_indices": [], "passthrough_row_indices": all_zarr_rows,
                "eigenvalues": [], "k_stored": 0,
            })
            total_passthrough += len(all_zarr_rows)
            continue

        G = G[:, keep]
        zarr_rows = [r for r, k in zip(zarr_rows, keep) if k]
        flip_flags = [f for f, k in zip(flip_flags, keep) if k]

        # Apply ALID flip: count copies of ALID A1 (flip dosage where needed)
        flip_arr = np.array(flip_flags)
        G[:, flip_arr] = 2.0 - G[:, flip_arr]

        # Remove zero-variance columns (would produce NaN in corrcoef)
        stdev = G.std(axis=0)
        nonzero = stdev > 0
        G = G[:, nonzero]
        zarr_rows = [r for r, k in zip(zarr_rows, nonzero) if k]
        flip_flags = [f for f, k in zip(flip_flags, nonzero) if k]

        if len(zarr_rows) < MIN_VARIANTS:
            block_manifest.append({
                "block_id": bi, "chrom": chrom, "start": start, "stop": stop,
                "n_1kg_variants": len(block_bim), "n_zarr_variants": len(all_zarr_rows),
                "n_intersection_variants": 0, "n_passthrough_variants": len(all_zarr_rows),
                "zarr_row_indices": [], "passthrough_row_indices": all_zarr_rows,
                "eigenvalues": [], "k_stored": 0,
            })
            total_passthrough += len(all_zarr_rows)
            continue

        # Pearson correlation matrix and eigendecomposition
        R = np.corrcoef(G.T)
        eigenvalues, V = np.linalg.eigh(R.astype(np.float64))
        # eigh returns ascending order; reverse to descending
        eigenvalues = eigenvalues[::-1].copy()
        V = V[:, ::-1].copy()

        # Truncate to VAR_STORE_THRESH
        cumvar = np.cumsum(eigenvalues) / eigenvalues.sum()
        k_stored = int(np.searchsorted(cumvar, VAR_STORE_THRESH)) + 1
        k_stored = min(k_stored, len(eigenvalues))

        eigenvectors[f"b{bi}"] = V[:, :k_stored].astype(np.float32)

        passthrough_rows_final = [r for r in all_zarr_rows if r not in zarr_rows_set or r not in zarr_rows]
        # recalculate: passthrough = all zarr in block minus those in intersection (after MAF filter)
        intersection_set = set(zarr_rows)
        passthrough_rows_final = [r for r in all_zarr_rows if r not in intersection_set]

        block_manifest.append({
            "block_id": bi,
            "chrom": chrom,
            "start": start,
            "stop": stop,
            "n_1kg_variants": len(block_bim),
            "n_zarr_variants": len(all_zarr_rows),
            "n_intersection_variants": len(zarr_rows),
            "n_passthrough_variants": len(passthrough_rows_final),
            "zarr_row_indices": zarr_rows,
            "passthrough_row_indices": passthrough_rows_final,
            "eigenvalues": eigenvalues.tolist(),
            "k_stored": k_stored,
        })

        total_compressed += len(zarr_rows)
        total_passthrough += len(passthrough_rows_final)
        print(f"  Block {bi:3d}: m={len(zarr_rows):4d} intersection, "
              f"K={k_stored:4d} stored, {len(passthrough_rows_final):4d} passthrough", end="\r")

    print(f"\n\nTotal: {total_compressed:,} compressed, {total_passthrough:,} passthrough")

    print("Saving eigenvectors.npz ...")
    np.savez_compressed(OUT_DIR / "eigenvectors.npz", **eigenvectors)

    with open(OUT_DIR / "block_manifest.json", "w") as fh:
        json.dump(block_manifest, fh)

    print(f"Done. Written: {OUT_DIR}")


if __name__ == "__main__":
    main()
