"""Query module for SQLite-indexed BESD data."""

import sqlite3
import math
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np


def norm_cdf(z: float) -> float:
    """Approximate normal CDF using Abramowitz and Stegun formula.

    Accurate to about 0.00012.
    """
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = 1 if z >= 0 else -1
    z = abs(z) / math.sqrt(2)

    t = 1.0 / (1.0 + p * z)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-z * z)

    return 0.5 * (1.0 + sign * y)


class BESDQueryIndex:
    """Query BESD data from SQLite index database."""

    def __init__(self, db_path: str):
        """Initialize query engine with SQLite database.

        Args:
            db_path: Path to SQLite database file
        """
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")

        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row
        self._load_metadata()

    def _load_metadata(self) -> None:
        """Load metadata from database."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT key, value FROM besd_meta")
        self.metadata = {row['key']: row['value'] for row in cursor.fetchall()}

    def close(self) -> None:
        """Close database connection."""
        self.conn.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def query_snp_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        """Query SNPs in chromosome range (coordinates in kb).

        Args:
            chr_val: Chromosome
            start_kb: Start position in kb
            end_kb: End position in kb

        Returns:
            List of SNP records as dicts
        """
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT row_idx, chr, snp_id, genetic_dist, bp, a1, a2, freq
            FROM esi
            WHERE chr = ? AND bp >= ? AND bp <= ?
            ORDER BY bp
        """, (chr_val, start_bp, end_bp))

        return [dict(row) for row in cursor.fetchall()]

    def query_probe_range(self, chr_val: str, start_kb: float, end_kb: float) -> List[Dict]:
        """Query probes in chromosome range (coordinates in kb).

        Args:
            chr_val: Chromosome
            start_kb: Start position in kb
            end_kb: End position in kb

        Returns:
            List of probe records as dicts
        """
        start_bp = int(start_kb * 1000)
        end_bp = int(end_kb * 1000)

        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT row_idx, chr, probe_id, genetic_dist, probe_bp, gene, orientation
            FROM epi
            WHERE chr = ? AND probe_bp >= ? AND probe_bp <= ?
            ORDER BY probe_bp
        """, (chr_val, start_bp, end_bp))

        return [dict(row) for row in cursor.fetchall()]

    def get_probe_snps(self, probe_idx: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get SNP indices, betas, and SEs for a probe.

        Args:
            probe_idx: Probe row index

        Returns:
            Tuple of (snp_indices, betas, ses) as numpy arrays
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT snp_indices, betas, ses, snp_count
            FROM probe_data
            WHERE probe_idx = ?
        """, (probe_idx,))

        row = cursor.fetchone()
        if not row:
            return np.array([], dtype=np.int32), np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        snp_count = row['snp_count']
        if snp_count == 0:
            return np.array([], dtype=np.int32), np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        # Deserialize BLOBs
        snp_indices = np.frombuffer(row['snp_indices'], dtype=np.int32)
        betas = np.frombuffer(row['betas'], dtype=np.float32)
        ses = np.frombuffer(row['ses'], dtype=np.float32)

        return snp_indices, betas, ses

    def query_cis_window(
        self,
        snp_chr: str, snp_start_kb: float, snp_end_kb: float,
        probe_chr: str, probe_start_kb: float, probe_end_kb: float,
    ) -> List[Dict]:
        """Query cis-window: SNP region + probe region.

        Args:
            snp_chr: SNP chromosome
            snp_start_kb: SNP region start (kb)
            snp_end_kb: SNP region end (kb)
            probe_chr: Probe chromosome
            probe_start_kb: Probe region start (kb)
            probe_end_kb: Probe region end (kb)

        Returns:
            List of associations with metadata and statistics
        """
        # Get SNPs and probes in ranges
        snps = self.query_snp_range(snp_chr, snp_start_kb, snp_end_kb)
        probes = self.query_probe_range(probe_chr, probe_start_kb, probe_end_kb)

        # Build SNP index set for fast lookup
        snp_indices_set = {s['row_idx'] for s in snps}
        snp_by_idx = {s['row_idx']: s for s in snps}

        # Query associations
        associations = []
        for probe in probes:
            probe_idx = probe['row_idx']
            snp_indices, betas, ses = self.get_probe_snps(probe_idx)

            # Find associations that match SNP range
            for i, snp_idx in enumerate(snp_indices):
                if snp_idx in snp_indices_set:
                    snp = snp_by_idx[snp_idx]
                    beta = float(betas[i])
                    se = float(ses[i])

                    # Calculate p-value
                    if se > 0:
                        z_score = abs(beta / se)
                        pval = 2 * (1 - norm_cdf(z_score))
                    else:
                        pval = 1.0

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

    def query_by_probe_id(self, probe_id: str) -> List[Dict]:
        """Query all associations for a specific probe.

        Args:
            probe_id: Probe ID to query

        Returns:
            List of associations for the probe
        """
        cursor = self.conn.cursor()

        # Find probe
        cursor.execute("""
            SELECT row_idx, chr, probe_id, probe_bp, gene
            FROM epi
            WHERE probe_id = ?
        """, (probe_id,))

        probe_row = cursor.fetchone()
        if not probe_row:
            return []

        probe = dict(probe_row)
        probe_idx = probe['row_idx']

        # Get probe data
        snp_indices, betas, ses = self.get_probe_snps(probe_idx)
        if len(snp_indices) == 0:
            return []

        # Get SNP metadata - convert numpy int32 to Python int for SQLite
        snp_indices_list = [int(idx) for idx in snp_indices]
        placeholders = ','.join('?' * len(snp_indices_list))
        cursor.execute(f"""
            SELECT row_idx, snp_id, chr, bp, a1, a2
            FROM esi
            WHERE row_idx IN ({placeholders})
        """, snp_indices_list)

        snp_by_idx = {row['row_idx']: dict(row) for row in cursor.fetchall()}

        # Build results
        associations = []
        for i, snp_idx in enumerate(snp_indices):
            snp_idx_int = int(snp_idx)
            if snp_idx_int in snp_by_idx:
                snp = snp_by_idx[snp_idx_int]
                beta = float(betas[i])
                se = float(ses[i])

                if se > 0:
                    z_score = abs(beta / se)
                    pval = 2 * (1 - norm_cdf(z_score))
                else:
                    pval = 1.0

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

    def query_by_snp_id(self, snp_id: str) -> List[Dict]:
        """Query all associations for a specific SNP.

        Args:
            snp_id: SNP ID to query

        Returns:
            List of associations for the SNP
        """
        cursor = self.conn.cursor()

        # Find SNP
        cursor.execute("""
            SELECT row_idx, chr, snp_id, bp, a1, a2
            FROM esi
            WHERE snp_id = ?
        """, (snp_id,))

        snp_row = cursor.fetchone()
        if not snp_row:
            return []

        snp = dict(snp_row)
        target_snp_idx = snp['row_idx']

        # Search all probes for associations with this SNP
        associations = []
        cursor.execute("SELECT row_idx FROM epi")
        for probe_row in cursor.fetchall():
            probe_idx = probe_row['row_idx']
            snp_indices, betas, ses = self.get_probe_snps(probe_idx)

            # Find this SNP in the probe's data - convert to int for comparison
            target_idx_int = int(target_snp_idx)
            snp_indices_int = [int(idx) for idx in snp_indices]

            if target_idx_int in snp_indices_int:
                match_idx = snp_indices_int.index(target_idx_int)
                beta = float(betas[match_idx])
                se = float(ses[match_idx])

                if se > 0:
                    z_score = abs(beta / se)
                    pval = 2 * (1 - norm_cdf(z_score))
                else:
                    pval = 1.0

                # Get probe metadata
                cursor.execute("""
                    SELECT chr, probe_id, probe_bp, gene
                    FROM epi
                    WHERE row_idx = ?
                """, (probe_idx,))

                probe_data = dict(cursor.fetchone())

                associations.append({
                    'snp_id': snp['snp_id'],
                    'snp_chr': snp['chr'],
                    'snp_bp': snp['bp'],
                    'a1': snp['a1'],
                    'a2': snp['a2'],
                    'probe_id': probe_data['probe_id'],
                    'probe_chr': probe_data['chr'],
                    'probe_bp': probe_data['probe_bp'],
                    'gene': probe_data['gene'],
                    'beta': beta,
                    'se': se,
                    'pval': pval,
                })

        return associations
