#!/usr/bin/env python3
"""Download canonical TSS positions for all human HGNC genes from Ensembl BioMart.

Writes resources/ensembl_gene_tss_hg38.tsv with columns: gene_name, chr, bp.
One row per HGNC symbol; bp is the median TSS across all annotated transcripts
for that gene on its most common chromosome.

Usage:
    python3 scripts/fetch_ensembl_gene_tss.py [--output PATH] [--release N]
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from urllib.parse import urlencode
from urllib.request import urlopen

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

QUERY_XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="chromosome_name"
            value="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"/>
    <Attribute name="hgnc_symbol"/>
    <Attribute name="chromosome_name"/>
    <Attribute name="transcription_start_site"/>
  </Dataset>
</Query>"""


def fetch(timeout: int = 300) -> list[tuple[str, str, int]]:
    """Return list of (hgnc_symbol, chrom, tss) tuples from BioMart."""
    print("Querying Ensembl BioMart…", file=sys.stderr)
    data = urlencode({"query": QUERY_XML}).encode()
    with urlopen(BIOMART_URL, data=data, timeout=timeout) as resp:
        lines = resp.read().decode().splitlines()

    rows = []
    for line in lines[1:]:  # skip header
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        symbol, chrom, tss = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if not symbol or not chrom or not tss:
            continue
        try:
            rows.append((symbol, chrom, int(tss)))
        except ValueError:
            continue
    print(f"  {len(rows):,} transcript TSS records downloaded", file=sys.stderr)
    return rows


def canonical_tss(rows: list[tuple[str, str, int]]) -> list[tuple[str, str, int]]:
    """Collapse to one canonical TSS per HGNC symbol.

    Picks the most common chromosome for each symbol (handles PAR genes), then
    takes the median TSS across all transcripts on that chromosome.
    """
    by_gene: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for symbol, chrom, tss in rows:
        by_gene[symbol].append((chrom, tss))

    result = []
    for symbol, entries in by_gene.items():
        chrom_counts = Counter(c for c, _ in entries)
        best_chrom = chrom_counts.most_common(1)[0][0]
        positions = [bp for c, bp in entries if c == best_chrom]
        result.append((symbol, best_chrom, int(median(positions))))

    # Sort: numeric chromosomes (1-22) then X, Y, MT, then by bp
    def _key(item):
        sym, chrom, bp = item
        try:
            return (0, int(chrom), bp, sym)
        except ValueError:
            return (1, {'X': 23, 'Y': 24, 'MT': 25}.get(chrom, 99), bp, sym)

    return sorted(result, key=_key)


def write(records: list[tuple[str, str, int]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene_name', 'chr', 'bp'])
        for symbol, chrom, bp in records:
            writer.writerow([symbol, chrom, bp])
    print(f"Wrote {len(records):,} genes to {path}", file=sys.stderr)


def main() -> None:
    repo_root = Path(__file__).parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output', default=str(repo_root / 'resources' / 'ensembl_gene_tss_hg38.tsv'),
        help='Output TSV path (default: resources/ensembl_gene_tss_hg38.tsv)',
    )
    parser.add_argument(
        '--timeout', type=int, default=300,
        help='HTTP timeout in seconds (default: 300)',
    )
    args = parser.parse_args()

    raw = fetch(timeout=args.timeout)
    records = canonical_tss(raw)
    write(records, Path(args.output))


if __name__ == '__main__':
    main()
