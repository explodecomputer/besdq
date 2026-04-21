"""Command-line interface for BESD query tool."""

import argparse
import sys
from pathlib import Path
from typing import List, Dict

from .besd_reader import BESDQueryEngine
from .builder import BESDIndexBuilder


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

    return chr_part, start_bp, end_bp


def write_output(associations: List[Dict], output_path: str, pval_threshold: float = 0.05):
    """Write results in SMR-compatible format."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter by p-value threshold
    filtered = [a for a in associations if a['pval'] < pval_threshold]

    with open(f"{output_path}.txt", 'w') as f:
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
                assoc['probe_id'],
                str(assoc['probe_chr']),
                str(assoc['probe_bp']),
                assoc['gene'] or 'NA',
                f"{assoc['beta']:.6f}",
                f"{assoc['se']:.6f}",
                f"{assoc['pval']:.6e}",
            ]) + '\n')

    print(f"Results written to {output_path}.txt")
    print(f"  Associations (before filter): {len(associations)}")
    print(f"  Associations (after p-value filter < {pval_threshold}): {len(filtered)}")


def main():
    parser = argparse.ArgumentParser(
        description='Direct BESD query tool (SMR-compatible interface)'
    )
    parser.add_argument('--beqtl-summary', required=True,
                        help='Path to BESD files (without extension)')
    parser.add_argument('--index',
                        help='Create SQLite index database at specified path (e.g. data/westra.db)')
    parser.add_argument('--query', type=float, default=0.05,
                        help='P-value threshold for filtering results')
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
                        help='SNP region as chr:pos or chr:start-end (e.g. 1:1191870 or 1:100000-2000000)')
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
                        help='Probe region as chr:pos or chr:start-end (e.g. 1:1140818 or 1:1000000-2000000)')
    parser.add_argument('--out',
                        help='Output file prefix')

    args = parser.parse_args()

    # Handle indexing mode
    if args.index:
        try:
            builder = BESDIndexBuilder(args.index)
            builder.build(args.beqtl_summary, force=False)
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

    # Parse chrpos arguments if provided
    snp_chr = args.snp_chr
    snp_start_kb = None
    snp_end_kb = None
    probe_chr = args.probe_chr
    probe_start_kb = None
    probe_end_kb = None

    if args.snp_chrpos:
        if args.snp_chr or args.from_snp_kb or args.to_snp_kb or args.from_snp_bp or args.to_snp_bp:
            print("Error: Cannot use --snp-chrpos with --snp-chr, --from-snp-kb, --to-snp-kb, --from-snp-bp, or --to-snp-bp", file=sys.stderr)
            sys.exit(1)
        try:
            snp_chr, start_bp, end_bp = parse_chrpos(args.snp_chrpos)
            snp_start_kb = start_bp / 1000.0
            snp_end_kb = end_bp / 1000.0
        except ValueError as e:
            print(f"Error parsing --snp-chrpos: {e}", file=sys.stderr)
            sys.exit(1)

    if args.probe_chrpos:
        if args.probe_chr or args.from_probe_kb or args.to_probe_kb or args.from_probe_bp or args.to_probe_bp:
            print("Error: Cannot use --probe-chrpos with --probe-chr, --from-probe-kb, --to-probe-kb, --from-probe-bp, or --to-probe-bp", file=sys.stderr)
            sys.exit(1)
        try:
            probe_chr, start_bp, end_bp = parse_chrpos(args.probe_chrpos)
            probe_start_kb = start_bp / 1000.0
            probe_end_kb = end_bp / 1000.0
        except ValueError as e:
            print(f"Error parsing --probe-chrpos: {e}", file=sys.stderr)
            sys.exit(1)

    # Validate and convert SNP position arguments (if not using chrpos)
    if snp_start_kb is None:
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

    # Validate and convert probe position arguments (if not using chrpos)
    if probe_start_kb is None:
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

    try:
        print(f"Loading BESD data from {args.beqtl_summary}...")
        engine = BESDQueryEngine(args.beqtl_summary)

        print(f"SNPs loaded: {len(engine.snps)}")
        print(f"Probes loaded: {len(engine.probes)}")
        print(f"Format: SPARSE_FILE_TYPE_{engine.besd.format_type}")

        print(f"\nQuerying cis-window:")
        print(f"  SNP:   {snp_chr}:{snp_start_kb*1000:.0f}-{snp_end_kb*1000:.0f} bp")
        print(f"  Probe: {probe_chr}:{probe_start_kb*1000:.0f}-{probe_end_kb*1000:.0f} bp")

        associations = engine.query_cis_window(
            snp_chr=snp_chr,
            snp_start_kb=snp_start_kb,
            snp_end_kb=snp_end_kb,
            probe_chr=probe_chr,
            probe_start_kb=probe_start_kb,
            probe_end_kb=probe_end_kb,
        )

        print(f"\nFound {len(associations)} associations")

        if associations:
            write_output(associations, args.out, pval_threshold=args.query)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
