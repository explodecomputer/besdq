"""SQLite database builder for BESD data."""

import sqlite3
import struct
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np

from .besd_reader import IndexReader, BESDReader


class BESDIndexBuilder:
    """Build SQLite index from BESD files."""

    def __init__(self, db_path: str):
        """Initialize builder with database path."""
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

    def build(self, besd_prefix: str, force: bool = False) -> None:
        """Build index database from BESD files.

        Args:
            besd_prefix: Path to BESD files (without extension)
            force: Overwrite existing database if True
        """
        if self.db_path.exists() and not force:
            raise FileExistsError(f"Database {self.db_path} already exists. Use force=True to overwrite.")

        if self.db_path.exists():
            self.db_path.unlink()

        # Load BESD files
        esi_path = f"{besd_prefix}.esi"
        epi_path = f"{besd_prefix}.epi"
        besd_path = f"{besd_prefix}.besd"

        print(f"Loading BESD files from {besd_prefix}...")
        snps = IndexReader.read_esi(esi_path)
        probes = IndexReader.read_epi(epi_path)
        besd = BESDReader(besd_path, len(probes))

        print(f"Loaded {len(snps)} SNPs and {len(probes)} probes")
        print(f"BESD format: SPARSE_FILE_TYPE_{besd.format_type}")

        # Create database and schema
        print(f"Creating database at {self.db_path}...")
        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()

        self._create_schema(cursor)

        # Load metadata
        print("Writing metadata...")
        self._write_metadata(cursor, {
            'format_type': besd.format_type,
            'n_snps': str(len(snps)),
            'n_probes': str(len(probes)),
            'besd_path': besd_path,
            'esi_path': esi_path,
            'epi_path': epi_path,
        })

        # Load SNP index
        print("Writing SNP index...")
        self._write_snps(cursor, snps)

        # Load probe index
        print("Writing probe index...")
        self._write_probes(cursor, probes)

        # Load probe data (statistics)
        print("Writing probe data...")
        self._write_probe_data(cursor, besd, len(probes))

        # Create indices
        print("Creating indices...")
        cursor.execute("CREATE INDEX idx_esi_chr_bp ON esi(chr, bp)")
        cursor.execute("CREATE INDEX idx_esi_snp_id ON esi(snp_id)")
        cursor.execute("CREATE INDEX idx_epi_chr_bp ON epi(chr, probe_bp)")
        cursor.execute("CREATE INDEX idx_epi_probe_id ON epi(probe_id)")

        conn.commit()
        conn.close()

        print(f"Database created successfully at {self.db_path}")

    def _create_schema(self, cursor: sqlite3.Cursor) -> None:
        """Create database schema."""
        # Metadata table
        cursor.execute("""
            CREATE TABLE besd_meta (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            )
        """)

        # SNP index table
        cursor.execute("""
            CREATE TABLE esi (
                row_idx INTEGER PRIMARY KEY,
                chr TEXT NOT NULL,
                snp_id TEXT NOT NULL,
                genetic_dist REAL,
                bp INTEGER NOT NULL,
                a1 TEXT,
                a2 TEXT,
                freq REAL
            )
        """)

        # Probe index table
        cursor.execute("""
            CREATE TABLE epi (
                row_idx INTEGER PRIMARY KEY,
                chr TEXT NOT NULL,
                probe_id TEXT NOT NULL,
                genetic_dist REAL,
                probe_bp INTEGER NOT NULL,
                gene TEXT,
                orientation TEXT
            )
        """)

        # Probe data (statistics) table
        cursor.execute("""
            CREATE TABLE probe_data (
                probe_idx INTEGER PRIMARY KEY,
                snp_count INTEGER NOT NULL,
                snp_indices BLOB NOT NULL,
                betas BLOB NOT NULL,
                ses BLOB NOT NULL
            )
        """)

    def _write_metadata(self, cursor: sqlite3.Cursor, metadata: Dict[str, str]) -> None:
        """Write metadata to database."""
        for key, value in metadata.items():
            cursor.execute(
                "INSERT INTO besd_meta (key, value) VALUES (?, ?)",
                (key, value)
            )

    def _write_snps(self, cursor: sqlite3.Cursor, snps: List[Dict]) -> None:
        """Write SNP index to database."""
        for snp in snps:
            cursor.execute("""
                INSERT INTO esi (row_idx, chr, snp_id, genetic_dist, bp, a1, a2, freq)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                snp['row_idx'],
                snp['chr'],
                snp['snp_id'],
                snp['genetic_dist'],
                snp['bp'],
                snp['a1'],
                snp['a2'],
                snp['freq'],
            ))

    def _write_probes(self, cursor: sqlite3.Cursor, probes: List[Dict]) -> None:
        """Write probe index to database."""
        for probe in probes:
            cursor.execute("""
                INSERT INTO epi (row_idx, chr, probe_id, genetic_dist, probe_bp, gene, orientation)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                probe['row_idx'],
                probe['chr'],
                probe['probe_id'],
                probe['genetic_dist'],
                probe['probe_bp'],
                probe['gene'],
                probe['orientation'],
            ))

    def _write_probe_data(self, cursor: sqlite3.Cursor, besd: BESDReader, n_probes: int) -> None:
        """Write probe statistics data to database."""
        for probe_idx in range(n_probes):
            # Get associations for this probe
            assocs = besd.get_probe_associations(probe_idx)

            if not assocs:
                # Empty probe
                snp_indices = np.array([], dtype=np.int32)
                betas = np.array([], dtype=np.float32)
                ses = np.array([], dtype=np.float32)
            else:
                # Extract arrays
                snp_indices_list, betas_list, ses_list = zip(*assocs)
                snp_indices = np.array(snp_indices_list, dtype=np.int32)
                betas = np.array(betas_list, dtype=np.float32)
                ses = np.array(ses_list, dtype=np.float32)

            # Serialize as BLOBs
            cursor.execute("""
                INSERT INTO probe_data (probe_idx, snp_count, snp_indices, betas, ses)
                VALUES (?, ?, ?, ?, ?)
            """, (
                probe_idx,
                len(assocs),
                snp_indices.tobytes(),
                betas.tobytes(),
                ses.tobytes(),
            ))

            if (probe_idx + 1) % 1000 == 0:
                print(f"  Wrote {probe_idx + 1} / {n_probes} probes")
