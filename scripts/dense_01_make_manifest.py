#!/usr/bin/env python3
"""Issue 01: Generate manifest TSV for the dense Zarr batch builder.

Scans a directory of VCF files and extracts trait_id, file_path, trait_name,
and n from the ##SAMPLE header line. Writes manifest.tsv to the output path.

Usage:
    python dense_01_make_manifest.py data/vcf-ukb manifest.tsv
"""

import re
import sys
from pathlib import Path


def parse_sample_header(vcf_path: Path) -> dict:
    """Extract fields from the ##SAMPLE header line of a VCF."""
    with __import__("gzip").open(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if not line.startswith("##SAMPLE"):
                continue
            # ##SAMPLE=<ID=UKB-b:10001,TotalControls=456349,TotalCases=6661,StudyType=CaseControl>
            fields = dict(re.findall(r"(\w+)=([^,>]+)", line))
            n_controls = int(fields.get("TotalControls", 0))
            n_cases = int(fields.get("TotalCases", 0))
            study_type = fields.get("StudyType", "Continuous")
            n = n_controls + n_cases if study_type == "CaseControl" else n_controls
            return {"n": n, "study_type": study_type}
    raise ValueError(f"No ##SAMPLE header found in {vcf_path}")


def trait_id_from_path(vcf_path: Path) -> str:
    # ukb-b-10001_chr1.vcf.gz -> ukb-b-10001
    return vcf_path.name.split("_chr")[0]


def main(vcf_dir: str, out_path: str) -> None:
    vcf_dir = Path(vcf_dir)
    vcfs = sorted(vcf_dir.glob("*.vcf.gz"))
    if not vcfs:
        sys.exit(f"No .vcf.gz files found in {vcf_dir}")

    rows = []
    for vcf in vcfs:
        trait_id = trait_id_from_path(vcf)
        info = parse_sample_header(vcf)
        rows.append({
            "trait_id": trait_id,
            "file_path": str(vcf.resolve()),
            "trait_name": trait_id,
            "n": info["n"],
        })

    out = Path(out_path)
    with out.open("w") as fh:
        fh.write("trait_id\tfile_path\ttrait_name\tn\n")
        for row in rows:
            fh.write(f"{row['trait_id']}\t{row['file_path']}\t{row['trait_name']}\t{row['n']}\n")

    print(f"Written {len(rows)} traits to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <vcf_dir> <out_manifest.tsv>")
    main(sys.argv[1], sys.argv[2])
