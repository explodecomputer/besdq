#!/usr/bin/env python3
"""
Direct BESD query tool — queries .besd/.esi/.epi files without database.

Supports both SPARSE_FILE_TYPE_3 and SPARSE_FILE_TYPE_3F formats.

Mimics SMR command-line interface:
    python simple_query.py \
        --beqtl-summary <prefix> \
        --snp-chr <chr> \
        --from-snp-kb <kb> \
        --to-snp-kb <kb> \
        --probe-chr <chr> \
        --from-probe-kb <kb> \
        --to-probe-kb <kb> \
        --out <output-prefix>
"""

import argparse
import sys
import math
from pathlib import Path
import struct
from typing import List, Tuple, Dict


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


def norm_cdf(z: float) -> float:
    """Approximate normal CDF using Abramowitz and Stegun formula.

    Accurate to about 0.00012.
    """
    # Constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of z
    sign = 1 if z >= 0 else -1
    z = abs(z) / math.sqrt(2)

    # Calculate t
    t = 1.0 / (1.0 + p * z)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-z * z)

    return 0.5 * (1.0 + sign * y)


class IndexReader:
    """Read .esi or .epi index files."""

    @staticmethod
    def read_esi(esi_path: str) -> List[Dict]:
        """Read SNP index file."""
        snps = []
        with open(esi_path, 'r') as f:
            for row_idx, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 4:
                    continue

                try:
                    snps.append({
                        'row_idx': row_idx,
                        'chr': parts[0],
                        'snp_id': parts[1],
                        'genetic_dist': float(parts[2]) if parts[2] != 'NA' else None,
                        'bp': int(parts[3]),
                        'a1': parts[4] if len(parts) > 4 else None,
                        'a2': parts[5] if len(parts) > 5 else None,
                        'freq': float(parts[6]) if len(parts) > 6 and parts[6] != 'NA' else None,
                    })
                except (ValueError, IndexError):
                    continue

        return snps

    @staticmethod
    def read_epi(epi_path: str) -> List[Dict]:
        """Read probe index file."""
        probes = []
        with open(epi_path, 'r') as f:
            for row_idx, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 4:
                    continue

                try:
                    probes.append({
                        'row_idx': row_idx,
                        'chr': parts[0],
                        'probe_id': parts[1],
                        'genetic_dist': float(parts[2]) if len(parts) > 2 and parts[2] != 'NA' else None,
                        'probe_bp': int(parts[3]),
                        'gene': parts[4] if len(parts) > 4 else None,
                        'orientation': parts[5] if len(parts) > 5 else None,
                    })
                except (ValueError, IndexError):
                    continue

        return probes


class BESDReader:
    """Read SPARSE_FILE_TYPE_3 and SPARSE_FILE_TYPE_3F BESD format files."""

    MAGIC_SPARSE_3F = 0x40400000  # Old sparse format
    MAGIC_SPARSE_3 = 3             # New sparse format
    RESERVEDUNITS = 16

    def __init__(self, besd_path: str, n_probes: int):
        """Initialize BESD reader."""
        self.besd_path = besd_path
        self.n_probes = n_probes
        self.format_type = None
        self._cols = None
        self._rowid = None
        self._val = None
        self._val_num = None
        self._open()

    def _open(self):
        """Open and parse BESD file."""
        with open(self.besd_path, 'rb') as f:
            # Read magic
            magic_bytes = f.read(4)
            magic = struct.unpack('<I', magic_bytes)[0]

            if magic == self.MAGIC_SPARSE_3F:
                self.format_type = '3F'
                self._parse_format_3f(f)
            elif magic == self.MAGIC_SPARSE_3:
                self.format_type = '3'
                self._parse_format_3(f)
            else:
                raise ValueError(f"Unsupported BESD magic: 0x{magic:08x}")

    def _parse_format_3f(self, f):
        """Parse SPARSE_FILE_TYPE_3F format."""
        # Read valNum (number of associations)
        val_num_bytes = f.read(8)
        self._val_num = struct.unpack('<Q', val_num_bytes)[0]

        # Calculate structure sizes
        col_num = (self.n_probes << 1) + 1  # probNum * 2 + 1

        # Read columns (offsets)
        col_bytes = f.read(col_num * 8)
        self._cols = list(struct.unpack(f'<{col_num}Q', col_bytes))

        # Read row IDs (SNP indices)
        row_id_bytes = f.read(self._val_num * 4)
        self._rowid = list(struct.unpack(f'<{self._val_num}I', row_id_bytes))

        # Read values (betas only - SEs are not in the BESD file)
        val_bytes = f.read(self._val_num * 4)
        self._val = list(struct.unpack(f'<{self._val_num}f', val_bytes))

    def _parse_format_3(self, f):
        """Parse SPARSE_FILE_TYPE_3 format."""
        # Read reserved units (sample size, n_snps, n_probes, etc.)
        reserved_bytes = f.read((self.RESERVEDUNITS - 1) * 4)
        reserved = struct.unpack(f'<{self.RESERVEDUNITS - 1}I', reserved_bytes)

        # Read valNum (number of associations)
        val_num_bytes = f.read(8)
        self._val_num = struct.unpack('<Q', val_num_bytes)[0]

        # Calculate structure sizes
        col_num = (self.n_probes << 1) + 1  # probNum * 2 + 1

        # Read columns (offsets)
        col_bytes = f.read(col_num * 8)
        self._cols = list(struct.unpack(f'<{col_num}Q', col_bytes))

        # Read row IDs (SNP indices)
        row_id_bytes = f.read(self._val_num * 4)
        self._rowid = list(struct.unpack(f'<{self._val_num}I', row_id_bytes))

        # Read values (betas only - SEs are not in the BESD file)
        val_bytes = f.read(self._val_num * 4)
        self._val = list(struct.unpack(f'<{self._val_num}f', val_bytes))

    def get_probe_associations(self, probe_idx: int) -> List[Tuple[int, float, float]]:
        """Get all associations for a probe.

        Returns: list of (snp_idx, beta, se) tuples

        Structure:
        - _cols[2p] = start of betas for probe p
        - _cols[2p+1] = start of SEs for probe p (= end of betas)
        - _cols[2p+2] = end of SEs for probe p
        """
        if probe_idx >= self.n_probes:
            return []

        if self._cols is None or self._rowid is None or self._val is None:
            return []

        # Get column indices for this probe
        beta_start = self._cols[probe_idx << 1]          # Start of betas
        se_start = self._cols[(probe_idx << 1) + 1]      # Start of SEs (= end of betas)
        se_end = self._cols[(probe_idx << 1) + 2]        # End of SEs

        num_assocs = se_start - beta_start

        if num_assocs <= 0:
            return []

        # Extract associations
        associations = []
        for i in range(int(num_assocs)):
            beta_idx = int(beta_start + i)
            se_idx = int(se_start + i)

            if beta_idx < len(self._rowid) and beta_idx < len(self._val) and se_idx < len(self._val):
                snp_idx = self._rowid[beta_idx]
                beta = self._val[beta_idx]
                se = self._val[se_idx]
                associations.append((snp_idx, beta, se))

        return associations


class BESDQueryEngine:
    """Query BESD data directly without database."""

    def __init__(self, besd_prefix: str):
        """Initialize from BESD file prefix (without extension)."""
        self.prefix = Path(besd_prefix)

        # Load indices
        esi_path = f"{besd_prefix}.esi"
        epi_path = f"{besd_prefix}.epi"
        besd_path = f"{besd_prefix}.besd"

        self.snps = IndexReader.read_esi(esi_path)
        self.probes = IndexReader.read_epi(epi_path)
        self.besd = BESDReader(besd_path, len(self.probes))

        # Build lookup maps
        self.snp_by_idx = {s['row_idx']: s for s in self.snps}
        self.probe_by_idx = {p['row_idx']: p for p in self.probes}

    def query_snp_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        """Query SNPs in chromosome range (coordinates in kb)."""
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        result = []
        for snp in self.snps:
            if snp['chr'] == chr_val and start_bp <= snp['bp'] <= end_bp:
                result.append(snp)
        return sorted(result, key=lambda x: x['bp'])

    def query_probe_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        """Query probes in chromosome range (coordinates in kb)."""
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        result = []
        for probe in self.probes:
            if probe['chr'] == chr_val and start_bp <= probe['probe_bp'] <= end_bp:
                result.append(probe)
        return sorted(result, key=lambda x: x['probe_bp'])

    def query_cis_window(
        self,
        snp_chr: str, snp_start_kb: float, snp_end_kb: float,
        probe_chr: str, probe_start_kb: float, probe_end_kb: float,
    ) -> List[Dict]:
        """Query cis-window: return all SNP-probe associations in both regions."""

        snps = self.query_snp_range(snp_chr, snp_start_kb, snp_end_kb)
        probes = self.query_probe_range(probe_chr, probe_start_kb, probe_end_kb)

        # Build SNP index set for fast lookup
        snp_indices_set = {s['row_idx'] for s in snps}

        # Find all associations
        associations = []
        for probe in probes:
            probe_idx = probe['row_idx']
            assocs = self.besd.get_probe_associations(probe_idx)

            for snp_idx, beta, se in assocs:
                if snp_idx in snp_indices_set:
                    snp = self.snp_by_idx[snp_idx]
                    # Calculate p-value: two-tailed test
                    if se > 0:
                        z_score = abs(beta / se)
                        pval = 2 * (1 - norm_cdf(z_score))
                    else:
                        pval = 1.0  # Unable to calculate p-value without SE

                    associations.append({
                        'snp_id': snp['snp_id'],
                        'snp_chr': snp['chr'],
                        'snp_bp': snp['bp'],
                        'a1': snp['a1'],
                        'a2': snp['a2'],
                        'probe_id': probe['probe_id'],
                        'probe_chr': probe['chr'],
                        'probe_bp': probe['probe_bp'],
                        'gene': probe['gene'],
                        'beta': beta,
                        'se': se,
                        'pval': pval,
                    })

        return associations


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
    parser.add_argument('--out', required=True,
                        help='Output file prefix')

    args = parser.parse_args()

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
