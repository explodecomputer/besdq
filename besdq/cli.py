"""Command-line interface for BESD query tool."""

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from typing import List, Dict

from .besd_reader import BESDQueryEngine
from .builder import BESDIndexBuilder
from .sqlite_query import BESDQueryIndex



def parse_chrpos(chrpos_str: str) -> tuple[str, int, int]:
    """Parse chr:pos or chr:start-end format and return (chr, start_bp, end_bp).

    Supports:
    - Single position: 1:1191870
    - Range: 1:100000-2000000
    """
    if ':' not in chrpos_str:
        raise ValueError(f"Invalid chrpos format: {chrpos_str}. Expected format: chr:pos or chr:start-end")

    chr_part, pos_part = chrpos_str.split(':', 1)

    # Check if it's a range or single position
    if '-' in pos_part:
        # Range format: start-end
        try:
            start_bp, end_bp = map(int, pos_part.split('-'))
        except ValueError:
            raise ValueError(f"Invalid position format in {chrpos_str}. Expected: start-end (integers)")
    else:
        # Single position format: pos
        try:
            pos = int(pos_part)
            start_bp = pos
            end_bp = pos
        except ValueError:
            raise ValueError(f"Invalid position format in {chrpos_str}. Expected: integer position")

    if start_bp > end_bp:
        raise ValueError(f"Invalid range in {chrpos_str}. Start position must be <= end position")

    return chr_part, start_bp, end_bp


def write_output(associations: List[Dict], output_path: str, pval_threshold: float = 0.05):
    """Write results in SMR-compatible format."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter by p-value threshold
    filtered = [a for a in associations if a['pval'] < pval_threshold]

    with open(output_path, 'w') as f:
        # Header
        f.write('\t'.join([
            'SNP', 'SNP_Chr', 'SNP_bp', 'A1', 'A2',
            'Probe', 'Probe_Chr', 'Probe_bp', 'Gene',
            'Beta', 'SE', 'P_value'
        ]) + '\n')

        # Data
        for assoc in filtered:
            f.write('\t'.join([
                assoc['snp_id'],
                str(assoc['snp_chr']),
                str(assoc['snp_bp']),
                assoc['a1'] or 'NA',
                assoc['a2'] or 'NA',
                assoc['trait_id'],
                str(assoc['trait_chr']) if assoc['trait_chr'] is not None else 'NA',
                str(assoc['trait_bp']) if assoc['trait_bp'] is not None else 'NA',
                assoc['gene'] or 'NA',
                f"{assoc['beta']:.6f}",
                f"{assoc['se']:.6f}",
                f"{assoc['pval']:.6e}",
            ]) + '\n')

    print(f"Results written to {output_path}")
    print(f"  Associations (before filter): {len(associations)}")
    print(f"  Associations (after p-value filter < {pval_threshold}): {len(filtered)}")


def _distribution(values):
    """Return distribution stats dict for a list of numbers (Nones skipped).

    Returns None if there are no non-None values.
    """
    vs = [v for v in values if v is not None]
    if not vs:
        return None
    n = len(vs)
    vs_sorted = sorted(vs)

    def pct(p):
        idx = (n - 1) * p / 100.0
        lo = int(idx)
        hi = min(lo + 1, n - 1)
        return vs_sorted[lo] + (vs_sorted[hi] - vs_sorted[lo]) * (idx - lo)

    return {
        'n': n,
        'min': vs_sorted[0],
        'p25': pct(25),
        'p50': pct(50),
        'p75': pct(75),
        'max': vs_sorted[-1],
        'mean': sum(vs) / n,
    }


_DIST_COL_W = 12  # column width for distribution stats


def _fmt_dist_row(label: str, dist) -> str:
    """Format one distribution row for the human-readable table."""
    w = _DIST_COL_W
    if dist is None:
        dash = f"{'—':>{w}}"
        return f"{label:<22}" + dash * 6
    return (
        f"{label:<22}"
        f"{int(round(dist['min'])):>{w},}"
        f"{int(round(dist['p25'])):>{w},}"
        f"{int(round(dist['p50'])):>{w},}"
        f"{int(round(dist['p75'])):>{w},}"
        f"{int(round(dist['max'])):>{w},}"
        f"{dist['mean']:>{w}.1f}"
    )


def _show_info_db(db_path: str, as_json: bool = False) -> None:
    """Print summary of a BESDQ SQLite index to stdout."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    # Schema version check: require n_cis column
    cols = {row[1] for row in conn.execute("PRAGMA table_info(epi)")}
    if 'n_cis' not in cols:
        conn.close()
        print(
            "Error: database schema is outdated — please rebuild with import-gwas-ssf",
            file=sys.stderr,
        )
        sys.exit(1)

    meta = {row['key']: row['value'] for row in conn.execute("SELECT key, value FROM besd_meta")}

    traits = [dict(row) for row in conn.execute(
        "SELECT trait_id, trait_name, gene, context, trait_chr, trait_bp, "
        "n_source_snps, n_cis, n_sig_trans_peaks, n_sig_trans_peaks_approx, n_sug_trans "
        "FROM epi ORDER BY row_idx"
    )]
    counts = {row['probe_idx']: row['snp_count']
              for row in conn.execute("SELECT probe_idx, snp_count FROM probe_data")}
    conn.close()

    total_assocs = sum(counts.values())
    n_traits = len(traits)
    study_meta = {}
    if 'study_metadata' in meta:
        try:
            study_meta = json.loads(meta['study_metadata'])
        except (json.JSONDecodeError, TypeError):
            pass

    # Compute distributions
    dist_stored      = _distribution([counts.get(i, 0) for i in range(n_traits)])
    dist_source_snps = _distribution([t['n_source_snps'] for t in traits])
    dist_cis         = _distribution([t['n_cis'] for t in traits])
    dist_sig_peaks   = _distribution([t['n_sig_trans_peaks'] for t in traits])
    dist_sug_trans   = _distribution([t['n_sug_trans'] for t in traits])

    if as_json:
        out = {
            'path': db_path,
            'source': meta.get('source', 'unknown'),
            'n_traits': n_traits,
            'n_snps': int(meta.get('n_snps', 0)),
            'total_associations': total_assocs,
            'build_params': {
                'cis_radius': int(meta['cis_radius']) if 'cis_radius' in meta else None,
                'sig_threshold': float(meta['sig_threshold']) if 'sig_threshold' in meta else None,
                'sug_threshold': float(meta['sug_threshold']) if 'sug_threshold' in meta else None,
                'sig_radius': int(meta['sig_radius']) if 'sig_radius' in meta else None,
            },
            'study_metadata': study_meta,
            'distributions': {
                'snp_count_stored': dist_stored,
                'n_source_snps': dist_source_snps,
                'n_cis': dist_cis,
                'n_sig_trans_peaks': dist_sig_peaks,
                'n_sug_trans': dist_sug_trans,
            },
        }
        print(json.dumps(out, indent=2))
        return

    # Human-readable output
    print(f"\nDataset:   {db_path}")
    print(f"Source:    {meta.get('source', 'unknown')}")
    print(f"Traits:    {n_traits:,}")
    print(f"SNPs:      {int(meta.get('n_snps', 0)):,}")
    print(f"Assocs:    {total_assocs:,}")

    if any(k in meta for k in ('cis_radius', 'sig_threshold', 'sug_threshold', 'sig_radius')):
        print("\nBuild parameters:")
        if 'cis_radius' in meta:
            print(f"  cis_radius:     {int(meta['cis_radius']):,} bp")
        if 'sig_threshold' in meta:
            print(f"  sig_threshold:  {meta['sig_threshold']}")
        if 'sug_threshold' in meta:
            print(f"  sug_threshold:  {meta['sug_threshold']}")
        if 'sig_radius' in meta:
            print(f"  sig_radius:     {int(meta['sig_radius']):,} bp")

    if study_meta:
        print("\nStudy metadata:")
        for key in ('genome_assembly', 'analysis_software', 'imputation_panel'):
            if key in study_meta:
                print(f"  {key}: {study_meta[key]}")
        if 'samples' in study_meta and study_meta['samples']:
            s = study_meta['samples'][0]
            if 'sample_size' in s:
                print(f"  sample_size: {s['sample_size']}")
            if 'sample_ancestry' in s:
                print(f"  ancestry: {', '.join(s['sample_ancestry'])}")

    # Distribution table
    w = _DIST_COL_W
    hdr = f"\n{'':22}{'min':>{w}}{'p25':>{w}}{'p50':>{w}}{'p75':>{w}}{'max':>{w}}{'mean':>{w}}"
    sep = "-" * (22 + w * 6)
    print(hdr)
    print(sep)
    print(_fmt_dist_row("snp_count (stored)", dist_stored))
    if dist_source_snps is not None:
        print(_fmt_dist_row("n_source_snps", dist_source_snps))
    if dist_cis is not None:
        print(_fmt_dist_row("n_cis", dist_cis))
    if dist_sig_peaks is not None:
        print(_fmt_dist_row("n_sig_trans_peaks", dist_sig_peaks))
    if dist_sug_trans is not None:
        print(_fmt_dist_row("n_sug_trans", dist_sug_trans))
    print()


def _show_info_besd(besd_path: str, as_json: bool = False) -> None:
    """Print summary of a legacy BESD binary dataset to stdout."""
    from .besd_reader import BESDQueryEngine
    engine = BESDQueryEngine(besd_path)

    n_probes = len(engine.probes)
    n_snps = len(engine.snps)

    # Total associations: _val_num counts betas+SEs together, so divide by 2
    total_assocs = (engine.besd._val_num // 2) if engine.besd._val_num else 0

    # Compute per-probe SNP counts from already-loaded cols array (O(1) per probe)
    snp_counts = []
    for probe in engine.probes:
        idx = probe['row_idx']
        if engine.besd._cols and (idx * 2 + 1) < len(engine.besd._cols):
            snp_counts.append(int(engine.besd._cols[idx * 2 + 1] - engine.besd._cols[idx * 2]))
        else:
            snp_counts.append(0)

    dist_stored = _distribution(snp_counts)

    if as_json:
        print(json.dumps({
            'path': besd_path,
            'source': 'besd',
            'n_probes': n_probes,
            'n_snps': n_snps,
            'total_associations': total_assocs,
            'distributions': {
                'snp_count_stored': dist_stored,
            },
        }, indent=2))
        return

    print(f"\nDataset:  {besd_path}")
    print(f"Source:   besd")
    print(f"Probes:   {n_probes:,}")
    print(f"SNPs:     {n_snps:,}")
    print(f"Assocs:   {total_assocs:,}")

    w = _DIST_COL_W
    hdr = f"\n{'':22}{'min':>{w}}{'p25':>{w}}{'p50':>{w}}{'p75':>{w}}{'max':>{w}}{'mean':>{w}}"
    print(hdr)
    print("-" * (22 + w * 6))
    print(_fmt_dist_row("snp_count (stored)", dist_stored))
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Direct BESD query tool (SMR-compatible interface)'
    )

    # Input source: either BESD files or SQLite index
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--besd',
                        help='Path to BESD files (without extension) - reads from binary files')
    input_group.add_argument('--db',
                        help='Path to SQLite index database file - reads from database')

    parser.add_argument('--build-db', dest='build_db',
                        help='Create SQLite index database at specified path (e.g. data/westra.db). Requires --besd.')
    parser.add_argument('--n-mode', choices=['scalar', 'vector'], default='scalar',
                        help='Sample-size storage mode: scalar (one n per probe) or vector (one n per SNP-probe pair). Default: scalar.')
    parser.add_argument('--sample-size', type=int, default=None,
                        help='Scalar sample size N used for --n-mode scalar. If omitted and ESI has allele frequencies, n is computed from the data.')
    parser.add_argument('--trait-variance', default=None, dest='trait_variance',
                        help='Two-column TSV (probe_id  var_y) supplying per-probe trait variance. Probes absent from the file default to var_y=1.0.')
    parser.add_argument('--query', type=float, default=1.0,
                        help='P-value threshold for filtering results (default: 1.0, i.e. no filter)')
    parser.add_argument('--original-scale', action='store_true', dest='original_scale',
                        help='Return beta and SE in original study units instead of SD units (default: SD units)')
    parser.add_argument('--snp-chr',
                        help='SNP chromosome')
    parser.add_argument('--from-snp-kb', type=float,
                        help='SNP region start (kb)')
    parser.add_argument('--to-snp-kb', type=float,
                        help='SNP region end (kb)')
    parser.add_argument('--from-snp-bp', type=int,
                        help='SNP region start (bp)')
    parser.add_argument('--to-snp-bp', type=int,
                        help='SNP region end (bp)')
    parser.add_argument('--snp-chrpos',
                        help='SNP region as chr:pos or chr:start-end (e.g. 1:1191870 or 1:100000-2000000). Comma-separated for multiple ranges.')
    parser.add_argument('--probe-chr',
                        help='Probe chromosome')
    parser.add_argument('--from-probe-kb', type=float,
                        help='Probe region start (kb)')
    parser.add_argument('--to-probe-kb', type=float,
                        help='Probe region end (kb)')
    parser.add_argument('--from-probe-bp', type=int,
                        help='Probe region start (bp)')
    parser.add_argument('--to-probe-bp', type=int,
                        help='Probe region end (bp)')
    parser.add_argument('--probe-chrpos',
                        help='Probe region as chr:pos or chr:start-end (e.g. 1:1140818 or 1:1000000-2000000). Comma-separated for multiple ranges.')
    parser.add_argument('--snp',
                        help='SNP ID(s) to query (comma-separated for multiple, e.g. rs123,rs456)')
    parser.add_argument('--probe',
                        help='Probe ID(s) to query (comma-separated for multiple)')
    parser.add_argument('--gene',
                        help='Gene name(s) to query (comma-separated for multiple)')
    parser.add_argument('--out',
                        help='Output file prefix')
    parser.add_argument('--info', action='store_true',
                        help='Print dataset summary (traits, SNPs, tier breakdown). Use with --db or --besd.')
    parser.add_argument('--json', action='store_true', dest='as_json',
                        help='Emit --info output as JSON instead of a human-readable table.')

    args = parser.parse_args()

    # Handle --info mode
    if args.info:
        if args.db:
            _show_info_db(args.db, as_json=args.as_json)
        else:
            _show_info_besd(args.besd, as_json=args.as_json)
        return

    # Handle indexing mode
    if args.build_db:
        if not args.besd:
            print("Error: --besd is required when using --build-db", file=sys.stderr)
            sys.exit(1)
        try:
            print(f"Building SQLite index from BESD files: {args.besd}")
            builder = BESDIndexBuilder(args.build_db)
            builder.build(
                args.besd,
                force=False,
                n_mode=args.n_mode,
                sample_size=args.sample_size,
                trait_variance_path=args.trait_variance,
            )
            print(f"Index created: {args.build_db}")
            return
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # Query mode requires --out
    if not args.out:
        print("Error: --out is required for query mode", file=sys.stderr)
        sys.exit(1)

    # Determine data source
    if args.db:
        print(f"Using SQLite index: {args.db}")
        query_engine = BESDQueryIndex(args.db, original_scale=args.original_scale)
        is_index = True
    else:
        print(f"Using BESD files: {args.besd}")
        query_engine = BESDQueryEngine(args.besd)
        is_index = False

    # Determine query type and validate arguments
    has_snp_query = bool(
        args.snp or args.snp_chrpos or args.snp_chr
        or args.from_snp_kb is not None or args.to_snp_kb is not None
        or args.from_snp_bp is not None or args.to_snp_bp is not None
    )
    has_probe_query = bool(
        args.probe or args.probe_chrpos or args.probe_chr or args.gene
        or args.from_probe_kb is not None or args.to_probe_kb is not None
        or args.from_probe_bp is not None or args.to_probe_bp is not None
    )
    id_query_count = sum(bool(x) for x in (args.snp, args.probe, args.gene))
    has_region_query_args = bool(
        args.snp_chrpos or args.probe_chrpos
        or args.snp_chr or args.probe_chr
        or args.from_snp_kb is not None or args.to_snp_kb is not None
        or args.from_snp_bp is not None or args.to_snp_bp is not None
        or args.from_probe_kb is not None or args.to_probe_kb is not None
        or args.from_probe_bp is not None or args.to_probe_bp is not None
    )

    if id_query_count > 1:
        print("Error: --snp, --probe, and --gene are mutually exclusive", file=sys.stderr)
        sys.exit(1)
    if id_query_count == 1 and has_region_query_args:
        print("Error: Cannot combine --snp/--probe/--gene with region query options", file=sys.stderr)
        sys.exit(1)
    if id_query_count == 0 and (not has_snp_query or not has_probe_query):
        print("Error: Cis-window queries require both SNP and probe query arguments", file=sys.stderr)
        sys.exit(1)

    if args.snp:
        # Query by SNP ID(s)
        snp_ids = [s.strip() for s in args.snp.split(',')]
        print(f"Querying {len(snp_ids)} SNP(s)...")
        associations = []
        for snp_id in snp_ids:
            assocs = query_engine.query_by_snp_id(snp_id)
            associations.extend(assocs)
            print(f"  {snp_id}: {len(assocs)} associations")

    elif args.probe:
        # Query by probe ID(s)
        probe_ids = [p.strip() for p in args.probe.split(',')]
        print(f"Querying {len(probe_ids)} probe(s)...")
        associations = []
        for probe_id in probe_ids:
            assocs = query_engine.query_by_probe_id(probe_id)
            associations.extend(assocs)
            print(f"  {probe_id}: {len(assocs)} associations")

    elif args.gene:
        # Query by gene name(s)
        genes = [g.strip() for g in args.gene.split(',')]
        print(f"Querying {len(genes)} gene(s)...")
        associations = []
        for gene_name in genes:
            assocs = query_engine.query_by_gene(gene_name)
            associations.extend(assocs)
            print(f"  {gene_name}: {len(assocs)} associations")

    else:
        # Cis-window query (original logic)
        # Handle comma-separated chrpos lists
        snp_ranges = []
        probe_ranges = []

        if args.snp_chrpos:
            if args.snp_chr or args.from_snp_kb or args.to_snp_kb or args.from_snp_bp or args.to_snp_bp:
                print("Error: Cannot use --snp-chrpos with --snp-chr, --from-snp-kb, --to-snp-kb, --from-snp-bp, or --to-snp-bp", file=sys.stderr)
                sys.exit(1)
            # Split by comma and parse each range
            try:
                for chrpos in args.snp_chrpos.split(','):
                    chrpos = chrpos.strip()
                    snp_chr, start_bp, end_bp = parse_chrpos(chrpos)
                    snp_ranges.append((snp_chr, start_bp / 1000.0, end_bp / 1000.0))
            except ValueError as e:
                print(f"Error parsing --snp-chrpos: {e}", file=sys.stderr)
                sys.exit(1)

        if args.probe_chrpos:
            if args.probe_chr or args.from_probe_kb or args.to_probe_kb or args.from_probe_bp or args.to_probe_bp:
                print("Error: Cannot use --probe-chrpos with --probe-chr, --from-probe-kb, --to-probe-kb, --from-probe-bp, or --to-probe-bp", file=sys.stderr)
                sys.exit(1)
            # Split by comma and parse each range
            try:
                for chrpos in args.probe_chrpos.split(','):
                    chrpos = chrpos.strip()
                    probe_chr, start_bp, end_bp = parse_chrpos(chrpos)
                    probe_ranges.append((probe_chr, start_bp / 1000.0, end_bp / 1000.0))
            except ValueError as e:
                print(f"Error parsing --probe-chrpos: {e}", file=sys.stderr)
                sys.exit(1)

        # If chrpos was used, combine results from all ranges
        if snp_ranges or probe_ranges:
            if not snp_ranges:
                print("Error: Must specify both --snp-chrpos and --probe-chrpos together", file=sys.stderr)
                sys.exit(1)
            if not probe_ranges:
                print("Error: Must specify both --snp-chrpos and --probe-chrpos together", file=sys.stderr)
                sys.exit(1)

            associations = []
            for snp_chr, snp_start_kb, snp_end_kb in snp_ranges:
                for probe_chr, probe_start_kb, probe_end_kb in probe_ranges:
                    print(f"Querying SNP {snp_chr}:{snp_start_kb*1000:.0f}-{snp_end_kb*1000:.0f} bp, Probe {probe_chr}:{probe_start_kb*1000:.0f}-{probe_end_kb*1000:.0f} bp")
                    range_assocs = query_engine.query_cis_window(
                        snp_chr=snp_chr,
                        snp_start_kb=snp_start_kb,
                        snp_end_kb=snp_end_kb,
                        probe_chr=probe_chr,
                        probe_start_kb=probe_start_kb,
                        probe_end_kb=probe_end_kb,
                    )
                    associations.extend(range_assocs)
        else:
            # Standard range query (not using chrpos)
            snp_chr = args.snp_chr
            snp_start_kb = None
            snp_end_kb = None
            probe_chr = args.probe_chr
            probe_start_kb = None
            probe_end_kb = None

            if args.from_snp_kb is not None and args.from_snp_bp is not None:
                print("Error: Cannot specify both --from-snp-kb and --from-snp-bp", file=sys.stderr)
                sys.exit(1)
            if args.to_snp_kb is not None and args.to_snp_bp is not None:
                print("Error: Cannot specify both --to-snp-kb and --to-snp-bp", file=sys.stderr)
                sys.exit(1)

            if args.from_snp_bp is not None:
                snp_start_kb = args.from_snp_bp / 1000.0
            elif args.from_snp_kb is not None:
                snp_start_kb = args.from_snp_kb
            else:
                print("Error: Must specify either --snp-chrpos, --from-snp-kb, or --from-snp-bp", file=sys.stderr)
                sys.exit(1)

            if args.to_snp_bp is not None:
                snp_end_kb = args.to_snp_bp / 1000.0
            elif args.to_snp_kb is not None:
                snp_end_kb = args.to_snp_kb
            else:
                print("Error: Must specify either --snp-chrpos, --to-snp-kb, or --to-snp-bp", file=sys.stderr)
                sys.exit(1)

            # Validate and convert probe position arguments
            if args.from_probe_kb is not None and args.from_probe_bp is not None:
                print("Error: Cannot specify both --from-probe-kb and --from-probe-bp", file=sys.stderr)
                sys.exit(1)
            if args.to_probe_kb is not None and args.to_probe_bp is not None:
                print("Error: Cannot specify both --to-probe-kb and --to-probe-bp", file=sys.stderr)
                sys.exit(1)

            if args.from_probe_bp is not None:
                probe_start_kb = args.from_probe_bp / 1000.0
            elif args.from_probe_kb is not None:
                probe_start_kb = args.from_probe_kb
            else:
                print("Error: Must specify either --probe-chrpos, --from-probe-kb, or --from-probe-bp", file=sys.stderr)
                sys.exit(1)

            if args.to_probe_bp is not None:
                probe_end_kb = args.to_probe_bp / 1000.0
            elif args.to_probe_kb is not None:
                probe_end_kb = args.to_probe_kb
            else:
                print("Error: Must specify either --probe-chrpos, --to-probe-kb, or --to-probe-bp", file=sys.stderr)
                sys.exit(1)

            print(f"Querying cis-window:")
            print(f"  SNP:   {snp_chr}:{snp_start_kb*1000:.0f}-{snp_end_kb*1000:.0f} bp")
            print(f"  Probe: {probe_chr}:{probe_start_kb*1000:.0f}-{probe_end_kb*1000:.0f} bp")

            associations = query_engine.query_cis_window(
                snp_chr=snp_chr,
                snp_start_kb=snp_start_kb,
                snp_end_kb=snp_end_kb,
                probe_chr=probe_chr,
                probe_start_kb=probe_start_kb,
                probe_end_kb=probe_end_kb,
            )

    try:
        print(f"\nFound {len(associations)} associations")

        if associations:
            write_output(associations, args.out, pval_threshold=args.query)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Clean up database connection if using index
        if is_index:
            query_engine.close()



def import_gwas_ssf_main() -> None:
    """Entry point for import-gwas-ssf command."""
    parser = argparse.ArgumentParser(
        prog='import-gwas-ssf',
        description='Import GWAS-SSF summary statistics into a queryable SQLite index',
    )
    parser.add_argument(
        '--trait-annotation', required=True, dest='trait_annotation',
        help='Path to trait annotation TSV file (required columns: file_path, trait_id, trait_name)',
    )
    parser.add_argument(
        '--ld-reference', required=False, default=None, dest='ld_reference',
        help='Prefix for plink2-format LD reference panel (--pfile prefix). '
             'Optional: if omitted, significant trans peaks are counted using a distance-based approximation.',
    )
    parser.add_argument(
        '--output', default=None,
        help='Output SQLite database path (default: <annotation_tsv_stem>.db)',
    )
    parser.add_argument(
        '--workers', type=int, default=1,
        help='Number of parallel workers for Pass 1 (default: 1)',
    )
    parser.add_argument(
        '--sig-threshold', type=float, default=5e-8, dest='sig_threshold',
        help='Genome-wide significance threshold (default: 5e-8)',
    )
    parser.add_argument(
        '--sug-threshold', type=float, default=1e-4, dest='sug_threshold',
        help='Suggestive significance threshold (default: 1e-4)',
    )
    parser.add_argument(
        '--cis-radius', type=int, default=1_000_000, dest='cis_radius',
        help='Cis window radius in bp (default: 1000000)',
    )
    parser.add_argument(
        '--sig-radius', type=int, default=500_000, dest='sig_radius',
        help='Significant trans window radius in bp (default: 500000)',
    )
    parser.add_argument(
        '--clump-r2', type=float, default=0.01, dest='clump_r2',
        help='LD clumping r² threshold (default: 0.01)',
    )
    parser.add_argument(
        '--clump-kb', type=int, default=10_000, dest='clump_kb',
        help='LD clumping window in kb (default: 10000)',
    )
    parser.add_argument(
        '--force-fallback', action='store_true', dest='force_fallback',
        help='Always use Python streaming reader instead of the pytabix+awk fast path',
    )

    args = parser.parse_args()

    # Route besdq.gwas_ssf_fast_reader INFO messages to stderr so the user can
    # see which read path is selected for each file.
    import logging as _logging
    _logging.basicConfig(
        level=_logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%H:%M:%S',
        stream=sys.stderr,
    )

    # Validate annotation file exists before starting processing
    annotation_path = Path(args.trait_annotation)
    if not annotation_path.exists():
        print(f"Error: annotation file not found: {args.trait_annotation}", file=sys.stderr)
        sys.exit(1)

    # Determine output path
    output_path = args.output or str(annotation_path.stem + '.db')

    try:
        from .annotation_reader import read_trait_annotation
        from .gwas_ssf_builder import GwasSsfIndexBuilder

        print(f"Reading trait annotations from: {args.trait_annotation}", file=sys.stderr)
        traits = read_trait_annotation(args.trait_annotation)
        print(f"  {len(traits)} traits found", file=sys.stderr)

        builder = GwasSsfIndexBuilder(output_path)
        builder.build(
            traits,
            workers=args.workers,
            cis_radius=args.cis_radius,
            sig_threshold=args.sig_threshold,
            sug_threshold=args.sug_threshold,
            plink2_pfile=args.ld_reference,
            sig_radius=args.sig_radius,
            clump_r2=args.clump_r2,
            clump_kb=args.clump_kb,
            force_fallback=args.force_fallback,
        )
        print(f"Index written to: {output_path}", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def discover_study_main() -> None:
    """Entry point for discover-study command."""
    parser = argparse.ArgumentParser(
        prog='besdq-discover-study',
        description='Query EBI GWAS Catalog for all studies under a PubMed ID',
    )
    parser.add_argument('pmid', help='PubMed ID (e.g. 38714679)')
    parser.add_argument(
        '--output', default=None,
        help='Write TSV to this file instead of stdout',
    )
    parser.add_argument(
        '--parse-gene', action='store_true', dest='parse_gene',
        help='Extract gene symbol from trait name and map to chr/bp',
    )
    parser.add_argument(
        '--gene-annotation', default=None, dest='gene_annotation',
        help='Path to gene annotation TSV (columns: gene_name, chr, bp). Required with --parse-gene.',
    )
    args = parser.parse_args()

    try:
        from .discover import discover_study, write_discovery_tsv
        rows = discover_study(
            args.pmid,
            parse_gene=args.parse_gene,
            gene_annotation_path=args.gene_annotation,
        )
        write_discovery_tsv(rows, output=args.output)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def build_stage1_main() -> None:
    """Entry point for build-stage1 command."""
    parser = argparse.ArgumentParser(
        prog='besdq-build-stage1',
        description='Filter GWAS-SSF source files and write intermediate TSV+YAML pairs',
    )
    parser.add_argument('discovery_tsv', help='Discovery TSV from discover-study')
    parser.add_argument('--workdir', required=True, help='Working directory for intermediates')
    parser.add_argument(
        '--parallel-downloads', type=int, default=1, dest='parallel_downloads',
        help='Number of concurrent download+filter workers (default: 1)',
    )
    parser.add_argument(
        '--sig-threshold', type=float, default=5e-8, dest='sig_threshold',
        help='Genome-wide significance threshold (default: 5e-8)',
    )
    parser.add_argument(
        '--sug-threshold', type=float, default=1e-4, dest='sug_threshold',
        help='Suggestive significance threshold (default: 1e-4)',
    )
    parser.add_argument(
        '--cis-radius', type=int, default=1_000_000, dest='cis_radius',
        help='Cis window radius in bp (default: 1000000)',
    )
    parser.add_argument(
        '--sig-radius', type=int, default=500_000, dest='sig_radius',
        help='Significant trans window radius in bp (default: 500000)',
    )
    parser.add_argument(
        '--ld-reference', default=None, dest='ld_reference',
        help='Prefix for plink2-format LD reference panel',
    )
    parser.add_argument(
        '--clump-r2', type=float, default=0.01, dest='clump_r2',
        help='LD clumping r² threshold (default: 0.01)',
    )
    parser.add_argument(
        '--clump-kb', type=int, default=10_000, dest='clump_kb',
        help='LD clumping window in kb (default: 10000)',
    )
    parser.add_argument(
        '--force-fallback', action='store_true', dest='force_fallback',
        help='Always use Python streaming reader',
    )
    args = parser.parse_args()

    try:
        from .stage1 import Stage1Builder
        builder = Stage1Builder()
        builder.build(
            discovery_tsv=args.discovery_tsv,
            workdir=args.workdir,
            parallel_downloads=args.parallel_downloads,
            cis_radius=args.cis_radius,
            sig_threshold=args.sig_threshold,
            sug_threshold=args.sug_threshold,
            plink2_pfile=args.ld_reference,
            sig_radius=args.sig_radius,
            clump_r2=args.clump_r2,
            clump_kb=args.clump_kb,
            force_fallback=args.force_fallback,
        )
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


def build_stage2_main() -> None:
    """Entry point for build-stage2 command."""
    parser = argparse.ArgumentParser(
        prog='besdq-build-stage2',
        description='Collate Stage 1 intermediate files into a queryable SQLite index',
    )
    parser.add_argument(
        'workdir_pmid',
        help='Path to the workdir/<pmid>/ directory containing intermediate files',
    )
    parser.add_argument('output_db', help='Output SQLite database path')
    args = parser.parse_args()

    try:
        from .gwas_ssf_builder import GwasSsfIndexBuilder
        builder = GwasSsfIndexBuilder(args.output_db)
        builder.build_from_intermediates(args.workdir_pmid)
        print(f"Index written to: {args.output_db}", file=sys.stderr)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
