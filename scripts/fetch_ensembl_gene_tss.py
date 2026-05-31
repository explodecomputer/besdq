#!/usr/bin/env python3
"""Download canonical TSS positions for all human HGNC genes and expand with aliases.

Writes resources/ensembl_gene_tss_hg38.tsv with columns:
  gene_name   – HGNC symbol (canonical) or alias / previous symbol
  chr         – chromosome (1-22, X, Y, MT)
  bp          – TSS position (GRCh38)
  canonical   – canonical HGNC symbol (empty for canonical rows, set for alias rows)

One canonical row per HGNC symbol, plus additional rows for each unambiguous
alias and previous symbol so that case-insensitive lookups on trait descriptions
resolve to the correct coordinates.

Usage:
    python3 scripts/fetch_ensembl_gene_tss.py [--output PATH] [--no-aliases]
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
HGNC_URL = (
    "https://www.genenames.org/cgi-bin/download/custom"
    "?col=gd_app_sym&col=gd_aliases&col=gd_prev_sym"
    "&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort"
    "&format=text&submit=submit"
)

BIOMART_QUERY = """<?xml version="1.0" encoding="UTF-8"?>
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


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------

def fetch_biomart(timeout: int = 300) -> list[tuple[str, str, int]]:
    """Return (hgnc_symbol, chrom, tss) tuples from Ensembl BioMart."""
    import requests as _req
    print("Fetching Ensembl BioMart TSS data…", file=sys.stderr)
    resp = _req.post(BIOMART_URL, data={"query": BIOMART_QUERY}, timeout=timeout)
    resp.raise_for_status()
    lines = resp.text.splitlines()
    rows = []
    for line in lines[1:]:
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        symbol, chrom, tss = parts[0].strip(), parts[1].strip(), parts[2].strip()
        if symbol and chrom and tss:
            try:
                rows.append((symbol, chrom, int(tss)))
            except ValueError:
                continue
    print(f"  {len(rows):,} transcript records", file=sys.stderr)
    return rows


def fetch_hgnc_aliases(timeout: int = 60) -> dict[str, dict]:
    """Return {canonical: {alias_symbols: [...], prev_symbols: [...]}} from HGNC."""
    import requests as _req
    print("Fetching HGNC alias data…", file=sys.stderr)
    resp = _req.get(HGNC_URL, timeout=timeout)
    resp.raise_for_status()
    lines = resp.text.splitlines()
    result: dict[str, dict] = {}
    reader = csv.DictReader(lines, delimiter='\t')
    for row in reader:
        canonical = row.get('Approved symbol', '').strip()
        if not canonical:
            continue
        aliases = [s.strip() for s in row.get('Alias symbols', '').split(',') if s.strip()]
        prev = [s.strip() for s in row.get('Previous symbols', '').split(',') if s.strip()]
        result[canonical] = {'alias': aliases, 'prev': prev}
    print(f"  {len(result):,} approved symbols with alias data", file=sys.stderr)
    return result


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def canonical_tss(rows: list[tuple[str, str, int]]) -> dict[str, tuple[str, int]]:
    """Collapse to one canonical TSS per HGNC symbol → {symbol: (chrom, bp)}."""
    by_gene: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for symbol, chrom, tss in rows:
        by_gene[symbol].append((chrom, tss))
    result = {}
    for symbol, entries in by_gene.items():
        best_chrom = Counter(c for c, _ in entries).most_common(1)[0][0]
        positions = [bp for c, bp in entries if c == best_chrom]
        result[symbol] = (best_chrom, int(median(positions)))
    return result


def expand_with_aliases(
    tss_map: dict[str, tuple[str, int]],
    hgnc: dict[str, dict],
) -> list[tuple[str, str, int, str]]:
    """Return (gene_name, chr, bp, canonical) rows for canonical + unambiguous aliases.

    Alias resolution priority:
      1. current alias_symbols  (preferred over prev_symbols when there's a conflict)
      2. prev_symbols

    Aliases that map to more than one canonical gene are skipped (ambiguous).
    Aliases that are already canonical gene names are skipped (would overwrite).
    """
    # Build alias → canonical mapping with conflict detection
    # alias_source: 'alias' beats 'prev' when the same alias appears for two genes
    alias_to_canonical: dict[str, tuple[str, str]] = {}  # alias_lower → (canonical, source)
    conflicts: set[str] = set()

    def _register(alias: str, canonical: str, source: str) -> None:
        key = alias.lower()
        if key in conflicts:
            return
        if key in alias_to_canonical:
            existing_canonical, existing_source = alias_to_canonical[key]
            if existing_canonical == canonical:
                return
            # Conflict: same alias → two different canonical genes
            # 'alias' source beats 'prev' source; same-priority conflict → drop
            if source == 'alias' and existing_source == 'prev':
                alias_to_canonical[key] = (canonical, source)
            elif source == 'prev' and existing_source == 'alias':
                pass  # keep existing
            else:
                conflicts.add(key)
                del alias_to_canonical[key]
        else:
            alias_to_canonical[key] = (canonical, source)

    for canonical, data in hgnc.items():
        if canonical not in tss_map:
            continue
        for alias in data['alias']:
            _register(alias, canonical, 'alias')
        for prev in data['prev']:
            _register(prev, canonical, 'prev')

    # Remove aliases that are canonical gene names (don't shadow them)
    canonical_lower = {s.lower() for s in tss_map}
    alias_to_canonical = {
        k: v for k, v in alias_to_canonical.items()
        if k not in canonical_lower
    }

    print(
        f"  {len(alias_to_canonical):,} unambiguous aliases "
        f"({len(conflicts):,} ambiguous skipped)",
        file=sys.stderr,
    )

    # Build output rows: canonical rows first, then alias rows
    def _chrom_key(chrom: str) -> tuple:
        try:
            return (0, int(chrom), 0)
        except ValueError:
            return (1, {'X': 23, 'Y': 24, 'MT': 25}.get(chrom, 99), 0)

    canonical_rows = [
        (symbol, chrom, bp, '')
        for symbol, (chrom, bp) in sorted(
            tss_map.items(),
            key=lambda x: (_chrom_key(x[1][0]), x[1][1], x[0]),
        )
    ]

    # Alias rows: sorted by alias name for determinism
    alias_rows = []
    for alias_lower, (canonical, _) in sorted(alias_to_canonical.items()):
        # Preserve original case from HGNC data by finding the alias in hgnc entries
        chrom, bp = tss_map[canonical]
        # Find original-case alias name from HGNC
        alias_original = alias_lower  # fallback
        for data in hgnc.values():
            for a in data['alias'] + data['prev']:
                if a.lower() == alias_lower:
                    alias_original = a
                    break
            else:
                continue
            break
        alias_rows.append((alias_original, chrom, bp, canonical))

    return canonical_rows + alias_rows


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def write(records: list[tuple[str, str, int, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene_name', 'chr', 'bp', 'canonical'])
        for name, chrom, bp, canonical in records:
            writer.writerow([name, chrom, bp, canonical])
    n_canonical = sum(1 for _, _, _, c in records if not c)
    n_alias = sum(1 for _, _, _, c in records if c)
    print(
        f"Wrote {len(records):,} rows to {path} "
        f"({n_canonical:,} canonical, {n_alias:,} aliases)",
        file=sys.stderr,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    repo_root = Path(__file__).parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output',
        default=str(repo_root / 'resources' / 'ensembl_gene_tss_hg38.tsv'),
        help='Output TSV path (default: resources/ensembl_gene_tss_hg38.tsv)',
    )
    parser.add_argument(
        '--no-aliases', action='store_true',
        help='Skip alias expansion (write canonical symbols only)',
    )
    parser.add_argument('--timeout', type=int, default=300)
    args = parser.parse_args()

    raw = fetch_biomart(timeout=args.timeout)
    tss_map = canonical_tss(raw)
    print(f"  {len(tss_map):,} canonical genes", file=sys.stderr)

    if args.no_aliases:
        records = [
            (sym, chrom, bp, '')
            for sym, (chrom, bp) in tss_map.items()
        ]
    else:
        hgnc = fetch_hgnc_aliases(timeout=args.timeout)
        records = expand_with_aliases(tss_map, hgnc)

    write(records, Path(args.output))


if __name__ == '__main__':
    main()
