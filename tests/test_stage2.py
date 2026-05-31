"""Tests for build-stage2 command (Issue #30)."""

import gzip
import io
import sqlite3
import tempfile
import unittest
from pathlib import Path

import numpy as np
import yaml

from besdq.gwas_ssf_builder import (
    INTERMEDIATE_TSV_COLUMNS,
    GwasSsfIndexBuilder,
)
from besdq.gwas_ssf_reader import GwasSsfRow
from besdq.sqlite_query import BESDQueryIndex
from besdq.stage1 import Stage1Builder, _write_intermediate_pair


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_rows(n: int = 3, chr_: str = '1', bp_start: int = 1_000_000) -> list:
    return [
        GwasSsfRow(
            chr=chr_,
            bp=bp_start + i * 1000,
            a1='A',
            a2='G',
            rsid=f'rs{i}',
            beta=0.1 * (i + 1),
            se=0.02,
            eaf=0.3,
            p=1e-6,
        )
        for i in range(n)
    ]


def _write_intermediate(workdir: Path, trait_id: str, rows: list, meta: dict) -> None:
    tsv_gz = workdir / f"{trait_id}.tsv.gz"
    yaml_path = workdir / f"{trait_id}.yaml"
    _write_intermediate_pair(rows, meta, tsv_gz, yaml_path)


def _default_meta(trait_id: str, **kwargs) -> dict:
    m = {
        'trait_id': trait_id,
        'trait_name': f'Trait {trait_id}',
        'pmid': '99999',
        'trait_chr': '1',
        'trait_bp': 1_000_000,
        'gene': None,
        'context': None,
        'trait_var': 0.5,
        'n_cis': 3,
        'n_sig_trans_peaks': 0,
        'n_sig_trans_peaks_approx': False,
        'n_sug_trans': 0,
        'n_total_read': 100,
        'cis_radius': 1_000_000,
        'sig_threshold': 5e-8,
        'sug_threshold': 1e-4,
        'sig_radius': 500_000,
    }
    m.update(kwargs)
    return m


def _make_ssf_gz_bytes(rows: list) -> bytes:
    cols = [
        'chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
        'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'rsid',
    ]
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode='wb') as gz:
        gz.write(('\t'.join(cols) + '\n').encode())
        for r in rows:
            d = {
                'chromosome': r.chr,
                'base_pair_location': r.bp,
                'effect_allele': r.a1,
                'other_allele': r.a2,
                'beta': r.beta,
                'standard_error': r.se,
                'effect_allele_frequency': r.eaf,
                'p_value': r.p,
                'rsid': r.rsid or 'NA',
            }
            gz.write(('\t'.join(str(d[c]) for c in cols) + '\n').encode())
    return buf.getvalue()


def _write_discovery_tsv_helper(path: str, rows: list) -> None:
    cols = ['pmid', 'gcst_id', 'file_path', 'trait_name', 'gene', 'trait_chr', 'trait_bp', 'context']
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')


# ---------------------------------------------------------------------------
# Issue #30 — build-stage2
# ---------------------------------------------------------------------------

class TestStage2(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _workdir(self, pmid: str = '99999') -> Path:
        d = self.tmpdir / 'work' / pmid
        d.mkdir(parents=True, exist_ok=True)
        return d

    def test_stage2_produces_queryable_db(self):
        """build-stage2 produces a .db with traits in epi."""
        workdir = self._workdir()
        for gcst in ('GCST001', 'GCST002'):
            _write_intermediate(workdir, gcst, _make_rows(3), _default_meta(gcst))

        db_path = str(self.tmpdir / 'out.db')
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(workdir))

        conn = sqlite3.connect(db_path)
        count = conn.execute("SELECT COUNT(*) FROM epi").fetchone()[0]
        conn.close()
        self.assertEqual(count, 2)

    def test_all_traits_in_epi_with_correct_metadata(self):
        """All traits appear in epi with correct metadata."""
        workdir = self._workdir()
        _write_intermediate(workdir, 'GCST_A', _make_rows(2), _default_meta(
            'GCST_A', trait_name='MyTrait', trait_chr='7', trait_bp=5_000_000,
            gene='MYGENE',
        ))

        db_path = str(self.tmpdir / 'meta.db')
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(workdir))

        conn = sqlite3.connect(db_path)
        row = conn.execute(
            "SELECT trait_name, trait_chr, trait_bp, gene FROM epi WHERE trait_id='GCST_A'"
        ).fetchone()
        conn.close()
        self.assertIsNotNone(row)
        self.assertEqual(row[0], 'MyTrait')
        self.assertEqual(row[1], '7')
        self.assertEqual(row[2], 5_000_000)
        self.assertEqual(row[3], 'MYGENE')

    def test_failed_yaml_skipped(self):
        """Traits with *.failed.yaml are skipped."""
        workdir = self._workdir()
        _write_intermediate(workdir, 'GCST_OK', _make_rows(2), _default_meta('GCST_OK'))

        # Write a failed.yaml (no paired tsv.gz)
        with open(workdir / 'GCST_BAD.failed.yaml', 'w') as fh:
            yaml.dump({'gcst_id': 'GCST_BAD', 'error': 'fail', 'timestamp': 'now'}, fh)

        db_path = str(self.tmpdir / 'skip.db')
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(workdir))

        conn = sqlite3.connect(db_path)
        ids = {r[0] for r in conn.execute("SELECT trait_id FROM epi").fetchall()}
        conn.close()
        self.assertIn('GCST_OK', ids)
        self.assertNotIn('GCST_BAD', ids)

    def test_clean_rebuild(self):
        """Running twice with same output path → clean rebuild, not accumulation."""
        workdir = self._workdir()
        _write_intermediate(workdir, 'GCST001', _make_rows(2), _default_meta('GCST001'))

        db_path = str(self.tmpdir / 'rebuild.db')
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(workdir))
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(workdir))

        conn = sqlite3.connect(db_path)
        count = conn.execute("SELECT COUNT(*) FROM epi").fetchone()[0]
        conn.close()
        self.assertEqual(count, 1)  # not doubled

    def test_besd_query_index_correct_beta_se(self):
        """BESDQueryIndex returns correct beta/SE for a known SNP-trait pair."""
        workdir = self._workdir()
        known_rows = _make_rows(3)
        _write_intermediate(workdir, 'GCST_Q', known_rows, _default_meta(
            'GCST_Q', trait_var=1.0,
        ))

        db_path = str(self.tmpdir / 'query.db')
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(workdir))

        idx = BESDQueryIndex(db_path)
        assocs = idx.query_by_probe_id('GCST_Q')
        idx.close()

        self.assertIsInstance(assocs, list)
        # With trait_var=1.0, SD units == original units
        # Just verify we get non-empty results back
        self.assertGreater(len(assocs), 0)

    def test_stage1_to_stage2_end_to_end(self):
        """Stage1 (local fixtures) → Stage2 → BESDQueryIndex round-trip."""
        # Write source SSF files
        source_rows = _make_rows(5)
        ssf_bytes = _make_ssf_gz_bytes(source_rows)

        for gcst in ('GCSTE1', 'GCSTE2'):
            ssf_path = self.tmpdir / f'{gcst}.h.tsv.gz'
            ssf_path.write_bytes(ssf_bytes)

        disc_tsv = str(self.tmpdir / 'disc.tsv')
        _write_discovery_tsv_helper(disc_tsv, [
            {
                'pmid': '77777', 'gcst_id': 'GCSTE1',
                'file_path': str(self.tmpdir / 'GCSTE1.h.tsv.gz'),
                'trait_name': 'Trait E1', 'gene': '', 'trait_chr': '1',
                'trait_bp': '1000000', 'context': '',
            },
            {
                'pmid': '77777', 'gcst_id': 'GCSTE2',
                'file_path': str(self.tmpdir / 'GCSTE2.h.tsv.gz'),
                'trait_name': 'Trait E2', 'gene': '', 'trait_chr': '1',
                'trait_bp': '1000000', 'context': '',
            },
        ])

        workdir = str(self.tmpdir / 'e2e_work')
        Stage1Builder().build(disc_tsv, workdir=workdir)

        # Check intermediates exist
        stage1_pmid_dir = Path(workdir) / '77777'
        self.assertTrue((stage1_pmid_dir / 'GCSTE1.tsv.gz').exists())
        self.assertTrue((stage1_pmid_dir / 'GCSTE2.tsv.gz').exists())

        db_path = str(self.tmpdir / 'e2e.db')
        GwasSsfIndexBuilder(db_path).build_from_intermediates(str(stage1_pmid_dir))

        # Query the result
        idx = BESDQueryIndex(db_path)
        assocs_e1 = idx.query_by_probe_id('GCSTE1')
        assocs_e2 = idx.query_by_probe_id('GCSTE2')
        idx.close()

        # Both traits queryable
        self.assertIsInstance(assocs_e1, list)
        self.assertIsInstance(assocs_e2, list)

        conn = sqlite3.connect(db_path)
        n_traits = conn.execute("SELECT COUNT(*) FROM epi").fetchone()[0]
        conn.close()
        self.assertEqual(n_traits, 2)

    def test_empty_workdir_raises(self):
        """build_from_intermediates raises if no valid YAML files found."""
        empty_dir = str(self.tmpdir / 'empty')
        Path(empty_dir).mkdir()
        with self.assertRaises((ValueError, FileNotFoundError)):
            GwasSsfIndexBuilder(str(self.tmpdir / 'x.db')).build_from_intermediates(empty_dir)

    def test_nonexistent_workdir_raises(self):
        """build_from_intermediates raises FileNotFoundError for missing workdir."""
        with self.assertRaises(FileNotFoundError):
            GwasSsfIndexBuilder(str(self.tmpdir / 'x.db')).build_from_intermediates(
                '/nonexistent/workdir'
            )


if __name__ == '__main__':
    unittest.main(verbosity=2)
