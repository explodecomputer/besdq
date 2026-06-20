#!/usr/bin/env python3
"""Issue 07: Download LDetect EUR LD block definitions and filter to chr1.

Writes data/ldetect_chr1.bed with columns: chr, start, stop.

Usage:
    python scripts/dense_07_get_ldetect.py
"""

import urllib.request
from pathlib import Path

URL = "https://raw.githubusercontent.com/MRCIEU/gwasglue/master/inst/extdata/ldetect/EUR.bed"
OUT = Path("data/ldetect_chr1.bed")


def main():
    print(f"Downloading LDetect EUR blocks from gwasglue...")
    with urllib.request.urlopen(URL) as resp:
        raw = resp.read().decode()

    lines = raw.strip().split("\n")
    header = lines[0]
    cols = header.split()
    print(f"  Header: {cols}")

    chr_col = next(i for i, c in enumerate(cols) if c.lower() in ("chr", "chrom"))
    start_col = next(i for i, c in enumerate(cols) if c.lower() == "start")
    stop_col = next(i for i, c in enumerate(cols) if c.lower() in ("stop", "end"))

    chr1_rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split()
        chrom = parts[chr_col].lstrip("chr")
        if chrom != "1":
            continue
        chr1_rows.append((chrom, int(parts[start_col]), int(parts[stop_col])))

    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("chr\tstart\tstop\n")
        for chrom, start, stop in chr1_rows:
            fh.write(f"{chrom}\t{start}\t{stop}\n")

    print(f"  Chr1 blocks: {len(chr1_rows)}")
    print(f"  Written: {OUT}")


if __name__ == "__main__":
    main()
