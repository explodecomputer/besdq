#!/usr/bin/env python3
"""Issue 02: Build union variant table from all VCFs in the manifest.

Single-pass scan: collects union of ALID keys (chr:pos:A1:A2, A1 alphabetically
first), accumulates EAF sum per variant to compute the average, and takes rsid
from the first file that sees each variant.

Writes a bgzipped TSV + tabix index inside the output .besdz directory.

Usage:
    python dense_02_build_variant_table.py manifest.tsv output.besdz
"""

import subprocess
import sys
from pathlib import Path

import numpy as np


BCFTOOLS = "/home/gh13047/miniforge3/envs/bcftools/bin/bcftools"
TABIX = "/home/gh13047/miniforge3/envs/bcftools/bin/tabix"
BGZIP = "/home/gh13047/miniforge3/envs/bcftools/bin/bgzip"


def alid(chrom: str, pos: int, ref: str, alt: str):
    """Return (A1, A2, flip) where A1 is alphabetically first allele.
    flip=True means the VCF effect allele (ALT) is A2 — caller must negate beta/Z.
    """
    if alt < ref:
        return alt, ref, False   # ALT is A1, effect allele is A1, no flip
    else:
        return ref, alt, True    # REF is A1, effect allele is A2, flip


def read_manifest(manifest_path: str) -> list[dict]:
    rows = []
    with open(manifest_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            parts = line.strip().split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def main(manifest_path: str, out_dir: str) -> None:
    manifest = read_manifest(manifest_path)
    out = Path(out_dir)
    out.mkdir(exist_ok=True)

    # key: (chrom, pos, a1, a2) -> [eaf_sum, count, rsid]
    variants: dict[tuple, list] = {}

    total_vcfs = len(manifest)
    for i, row in enumerate(manifest, 1):
        vcf_path = row["file_path"]
        print(f"[{i}/{total_vcfs}] {row['trait_id']}", end="\r", flush=True)

        cmd = [BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%AF]\n", vcf_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, pos_str, ref, alt, rsid = parts[0], parts[1], parts[2], parts[3], parts[4]
            if "," in alt:
                continue
            rsid = rsid if rsid and rsid != "." else "."
            try:
                af = float(parts[5])
            except ValueError:
                continue

            pos = int(pos_str)
            a1, a2, flip = alid(chrom, pos, ref, alt)
            key = (chrom, pos, a1, a2)

            # EAF of A1: if no flip, A1=ALT so eaf_a1=AF; if flip, A1=REF so eaf_a1=1-AF
            eaf_a1 = af if not flip else 1.0 - af

            if key not in variants:
                variants[key] = [eaf_a1, 1, rsid]
            else:
                variants[key][0] += eaf_a1
                variants[key][1] += 1
                if variants[key][2] == "." and rsid != ".":
                    variants[key][2] = rsid

        proc.wait()
        if proc.returncode != 0:
            print(f"\nWarning: bcftools failed (exit {proc.returncode}) for {vcf_path}", file=sys.stderr)

    print(f"\nUnion variant count: {len(variants):,}")

    # Sort by (chrom numerically then lexically, pos)
    def sort_key(k):
        chrom = k[0]
        try:
            return (0, int(chrom), k[1])
        except ValueError:
            return (1, chrom, k[1])

    sorted_keys = sorted(variants.keys(), key=sort_key)

    # Write TSV
    tsv_path = out / "variants.tsv"
    with open(tsv_path, "w") as fh:
        fh.write("chr\tpos\tA1\tA2\trsid\teaf\n")
        for key in sorted_keys:
            chrom, pos, a1, a2 = key
            eaf_sum, count, rsid = variants[key]
            eaf_avg = eaf_sum / count
            fh.write(f"{chrom}\t{pos}\t{a1}\t{a2}\t{rsid}\t{eaf_avg:.6f}\n")

    # bgzip and tabix
    gz_path = str(tsv_path) + ".gz"
    subprocess.run([BGZIP, "-f", str(tsv_path)], check=True)
    subprocess.run([TABIX, "-s1", "-b2", "-e2", "-S1", gz_path], check=True)
    print(f"Variant table written: {gz_path}")

    # Write row-index map (variant key -> row index) as numpy for fast lookup in builder
    keys_array = np.array([(c, p, a1, a2) for c, p, a1, a2 in sorted_keys], dtype=object)
    np.save(str(out / "variant_keys.npy"), keys_array)
    print(f"Variant key index written: {out / 'variant_keys.npy'}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <manifest.tsv> <output.besdz>")
    main(sys.argv[1], sys.argv[2])
