"""Tests for build-stage1 command (Issues #27, #28, #29)."""

import gzip
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml

from besdq.stage1 import (
    DiscoveryRow,
    Stage1Builder,
    _is_url,
    _write_failed_yaml,
    _write_intermediate_pair,
    read_discovery_tsv,
)
from besdq.gwas_ssf_builder import INTERMEDIATE_TSV_COLUMNS
from besdq.gwas_ssf_reader import GwasSsfRow


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_ssf_gz(rows: list) -> bytes:
    """Build a minimal GWAS-SSF gzip bytes for testing."""
    cols = [
        'chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
        'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'rsid',
    ]
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode='wb') as gz:
        gz.write(('\t'.join(cols) + '\n').encode())
        for r in rows:
            line = '\t'.join(str(r.get(c, 'NA')) for c in cols) + '\n'
            gz.write(line.encode())
    return buf.getvalue()


def _write_ssf_gz(rows: list, path: str) -> None:
    Path(path).write_bytes(_make_ssf_gz(rows))


def _make_cis_rows(n: int = 5, chr_: str = '1', bp_center: int = 1_000_000) -> list:
    """Make rows that will pass the cis filter."""
    return [
        {
            'chromosome': chr_,
            'base_pair_location': bp_center + i * 1000,
            'effect_allele': 'A',
            'other_allele': 'G',
            'beta': 0.1,
            'standard_error': 0.02,
            'effect_allele_frequency': 0.3,
            'p_value': 1e-6,
            'rsid': f'rs{i}',
        }
        for i in range(n)
    ]


def _write_discovery_tsv(path: str, rows: list) -> None:
    cols = ['pmid', 'gcst_id', 'file_path', 'trait_name', 'gene', 'trait_chr', 'trait_bp', 'context']
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')


# ---------------------------------------------------------------------------
# Issue #27 — local file tracer bullet
# ---------------------------------------------------------------------------

class TestStage1LocalFile(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _local_ssf(self, name: str, rows: list) -> str:
        path = str(Path(self.tmpdir) / name)
        _write_ssf_gz(rows, path)
        return path

    def test_local_path_writes_intermediate_pair(self):
        """Single-trait discovery TSV with local path → correct GCST*.tsv.gz + GCST*.yaml."""
        ssf_path = self._local_ssf('GCST001.h.tsv.gz', _make_cis_rows(5))
        discovery_tsv = str(Path(self.tmpdir) / 'discovery.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': 'GCST001',
            'file_path': ssf_path, 'trait_name': 'Test trait',
            'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])

        workdir = str(Path(self.tmpdir) / 'work')
        builder = Stage1Builder()
        counts = builder.build(discovery_tsv, workdir=workdir)

        self.assertEqual(counts['completed'], 1)
        self.assertEqual(counts['failed'], 0)

        out_dir = Path(workdir) / '99999'
        self.assertTrue((out_dir / 'GCST001.tsv.gz').exists())
        self.assertTrue((out_dir / 'GCST001.yaml').exists())

    def test_intermediate_tsv_contains_filtered_rows(self):
        """Intermediate TSV contains correct filtered rows."""
        ssf_path = self._local_ssf('GCST002.h.tsv.gz', _make_cis_rows(5))
        discovery_tsv = str(Path(self.tmpdir) / 'disc.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': 'GCST002',
            'file_path': ssf_path, 'trait_name': 'T2',
            'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])
        workdir = str(Path(self.tmpdir) / 'work2')
        Stage1Builder().build(discovery_tsv, workdir=workdir)

        tsv_gz = Path(workdir) / '99999' / 'GCST002.tsv.gz'
        with gzip.open(tsv_gz, 'rt') as fh:
            header = fh.readline().rstrip('\n').split('\t')
            data_lines = fh.readlines()

        self.assertEqual(header, INTERMEDIATE_TSV_COLUMNS)
        self.assertGreater(len(data_lines), 0)

    def test_yaml_sidecar_contains_required_fields(self):
        """YAML sidecar contains trait_chr, trait_bp, n_cis, n_sig_trans_peaks, n_sug_trans."""
        ssf_path = self._local_ssf('GCST003.h.tsv.gz', _make_cis_rows(5))
        discovery_tsv = str(Path(self.tmpdir) / 'disc3.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': 'GCST003',
            'file_path': ssf_path, 'trait_name': 'T3',
            'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])
        workdir = str(Path(self.tmpdir) / 'work3')
        Stage1Builder().build(discovery_tsv, workdir=workdir)

        yaml_path = Path(workdir) / '99999' / 'GCST003.yaml'
        with open(yaml_path) as fh:
            meta = yaml.safe_load(fh)

        for field in ('trait_chr', 'trait_bp', 'n_cis', 'n_sig_trans_peaks', 'n_sug_trans'):
            self.assertIn(field, meta, f"Missing field: {field}")

    def test_local_file_not_deleted(self):
        """Local file is NOT deleted after processing."""
        ssf_path = self._local_ssf('GCST004.h.tsv.gz', _make_cis_rows(5))
        discovery_tsv = str(Path(self.tmpdir) / 'disc4.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': 'GCST004',
            'file_path': ssf_path, 'trait_name': 'T4',
            'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])
        workdir = str(Path(self.tmpdir) / 'work4')
        Stage1Builder().build(discovery_tsv, workdir=workdir)

        self.assertTrue(Path(ssf_path).exists(), "Local source file was deleted — should not be")

    def test_atomic_write_no_partial_final_file(self):
        """Interrupted write leaves only .tmp, never a partial final file."""
        rows = [GwasSsfRow(chr='1', bp=1000, a1='A', a2='G', rsid='rs1',
                           beta=0.1, se=0.02, eaf=0.3, p=1e-6)]
        meta = {'trait_id': 'GCST_TMP', 'trait_chr': '1', 'trait_bp': 1000}
        tsv_gz = Path(self.tmpdir) / 'GCST_TMP.tsv.gz'
        yaml_path = Path(self.tmpdir) / 'GCST_TMP.yaml'

        # Simulate interrupted write: gzip write succeeds but yaml rename is killed
        original_replace = os.replace
        call_count = [0]

        def patched_replace(src, dst):
            call_count[0] += 1
            if call_count[0] == 2:  # second rename (yaml)
                raise OSError("Simulated kill")
            return original_replace(src, dst)

        with mock.patch('os.replace', side_effect=patched_replace):
            try:
                _write_intermediate_pair(rows, meta, tsv_gz, yaml_path)
            except OSError:
                pass

        # Final files: either both exist or neither (TSV.gz may have been renamed)
        # The YAML must not exist as a final file (was interrupted before its rename)
        self.assertFalse(yaml_path.exists(), "Partial YAML final file should not exist")

    def test_read_discovery_tsv_parses_correctly(self):
        """read_discovery_tsv returns correct DiscoveryRow objects."""
        tsv_path = str(Path(self.tmpdir) / 'test.tsv')
        _write_discovery_tsv(tsv_path, [
            {'pmid': '11111', 'gcst_id': 'GCST100', 'file_path': '/tmp/f.gz',
             'trait_name': 'Test', 'gene': 'IL10', 'trait_chr': '1',
             'trait_bp': '206941637', 'context': 'baseline'},
        ])
        rows = read_discovery_tsv(tsv_path)
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r.pmid, '11111')
        self.assertEqual(r.gcst_id, 'GCST100')
        self.assertEqual(r.trait_chr, '1')
        self.assertEqual(r.trait_bp, 206941637)
        self.assertEqual(r.gene, 'IL10')


# ---------------------------------------------------------------------------
# Issue #28 — URL download, delete-after-filter, --parallel-downloads
# ---------------------------------------------------------------------------

class TestStage1UrlDownload(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_url_detection(self):
        self.assertTrue(_is_url('https://example.com/file.gz'))
        self.assertTrue(_is_url('http://example.com/file.gz'))
        self.assertTrue(_is_url('ftp://example.com/file.gz'))
        self.assertFalse(_is_url('/local/path/file.gz'))
        self.assertFalse(_is_url('relative/path.gz'))

    def test_url_source_downloaded_filtered_and_deleted(self):
        """URL → downloaded, filtered, intermediate written, download deleted."""
        ssf_bytes = _make_ssf_gz(_make_cis_rows(5))
        discovery_tsv = str(Path(self.tmpdir) / 'disc.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': 'GCST010',
            'file_path': 'https://fake.ebi.ac.uk/GCST010.h.tsv.gz',
            'trait_name': 'T10', 'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])
        workdir = str(Path(self.tmpdir) / 'work')

        mock_resp = mock.MagicMock()
        mock_resp.status_code = 200
        mock_resp.iter_content.return_value = [ssf_bytes]

        with mock.patch('besdq.stage1._download_file') as mock_dl:
            def fake_download(url, dest):
                dest.write_bytes(ssf_bytes)
            mock_dl.side_effect = fake_download

            counts = Stage1Builder().build(discovery_tsv, workdir=workdir)

        self.assertEqual(counts['completed'], 1)
        # Intermediate pair exists
        self.assertTrue((Path(workdir) / '99999' / 'GCST010.tsv.gz').exists())
        self.assertTrue((Path(workdir) / '99999' / 'GCST010.yaml').exists())
        # Downloaded temp file deleted
        download_tmp = Path(workdir) / '99999' / 'GCST010.download.tmp'
        self.assertFalse(download_tmp.exists())

    def test_local_path_not_deleted_regression(self):
        """Local source locator is NOT deleted (regression from URL handling)."""
        ssf_path = str(Path(self.tmpdir) / 'GCST011.h.tsv.gz')
        _write_ssf_gz(_make_cis_rows(3), ssf_path)
        discovery_tsv = str(Path(self.tmpdir) / 'disc.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': 'GCST011',
            'file_path': ssf_path, 'trait_name': 'T11',
            'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])
        workdir = str(Path(self.tmpdir) / 'work')
        Stage1Builder().build(discovery_tsv, workdir=workdir)
        self.assertTrue(Path(ssf_path).exists())

    def test_parallel_downloads(self):
        """--parallel-downloads N runs multiple workers."""
        ssf_bytes = _make_ssf_gz(_make_cis_rows(3))
        discovery_rows = [
            {
                'pmid': '99999', 'gcst_id': f'GCST{i:03d}',
                'file_path': f'https://fake.ebi.ac.uk/GCST{i:03d}.h.tsv.gz',
                'trait_name': f'T{i}', 'gene': '', 'trait_chr': '1',
                'trait_bp': '1000000', 'context': '',
            }
            for i in range(4)
        ]
        discovery_tsv = str(Path(self.tmpdir) / 'disc.tsv')
        _write_discovery_tsv(discovery_tsv, discovery_rows)
        workdir = str(Path(self.tmpdir) / 'work')

        with mock.patch('besdq.stage1._download_file') as mock_dl:
            def fake_download(url, dest):
                dest.write_bytes(ssf_bytes)
            mock_dl.side_effect = fake_download

            counts = Stage1Builder().build(
                discovery_tsv, workdir=workdir, parallel_downloads=2
            )

        self.assertEqual(counts['completed'], 4)


# ---------------------------------------------------------------------------
# Issue #29 — Resumability
# ---------------------------------------------------------------------------

class TestStage1Resumability(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make_env(self, gcst_id: str = 'GCST001') -> tuple:
        ssf_path = str(Path(self.tmpdir) / f'{gcst_id}.h.tsv.gz')
        _write_ssf_gz(_make_cis_rows(3), ssf_path)
        discovery_tsv = str(Path(self.tmpdir) / 'discovery.tsv')
        _write_discovery_tsv(discovery_tsv, [{
            'pmid': '99999', 'gcst_id': gcst_id,
            'file_path': ssf_path, 'trait_name': 'T',
            'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
        }])
        workdir = str(Path(self.tmpdir) / 'work')
        return discovery_tsv, workdir

    def test_skip_completed_traits(self):
        """Re-running after complete run skips all traits."""
        disc, workdir = self._make_env()
        counts1 = Stage1Builder().build(disc, workdir=workdir)
        self.assertEqual(counts1['completed'], 1)

        counts2 = Stage1Builder().build(disc, workdir=workdir)
        self.assertEqual(counts2['skipped'], 1)
        self.assertEqual(counts2['completed'], 0)

    def test_retry_failed_traits(self):
        """Failed traits are retried on re-run; failed.yaml deleted on success."""
        disc, workdir = self._make_env('GCST002')
        work_pmid = Path(workdir) / '99999'
        work_pmid.mkdir(parents=True, exist_ok=True)

        # Pre-place a failed.yaml
        failed_path = work_pmid / 'GCST002.failed.yaml'
        _write_failed_yaml(failed_path, 'GCST002', 'Previous error')
        self.assertTrue(failed_path.exists())

        counts = Stage1Builder().build(disc, workdir=workdir)
        # Should retry and succeed
        self.assertEqual(counts['completed'], 1)
        # failed.yaml should be gone
        self.assertFalse(failed_path.exists())

    def test_tmp_file_treated_as_incomplete(self):
        """A .tmp file present but no final file → reprocessed on next run."""
        disc, workdir = self._make_env('GCST003')
        work_pmid = Path(workdir) / '99999'
        work_pmid.mkdir(parents=True, exist_ok=True)

        # Simulate a partial write: only .tmp exists
        (work_pmid / 'GCST003.tsv.gz.tmp').write_bytes(b'partial')

        counts = Stage1Builder().build(disc, workdir=workdir)
        self.assertEqual(counts['completed'], 1)
        self.assertTrue((work_pmid / 'GCST003.tsv.gz').exists())

    def test_failed_trait_does_not_abort_run(self):
        """A trait that fails → failed.yaml written, remaining traits still processed."""
        ssf_good = str(Path(self.tmpdir) / 'GCST_GOOD.h.tsv.gz')
        _write_ssf_gz(_make_cis_rows(3), ssf_good)

        disc_tsv = str(Path(self.tmpdir) / 'disc.tsv')
        _write_discovery_tsv(disc_tsv, [
            {
                'pmid': '99999', 'gcst_id': 'GCST_BAD',
                'file_path': '/nonexistent/path.tsv.gz',
                'trait_name': 'Bad', 'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
            },
            {
                'pmid': '99999', 'gcst_id': 'GCST_GOOD',
                'file_path': ssf_good,
                'trait_name': 'Good', 'gene': '', 'trait_chr': '1', 'trait_bp': '1000000', 'context': '',
            },
        ])
        workdir = str(Path(self.tmpdir) / 'work')
        counts = Stage1Builder().build(disc_tsv, workdir=workdir)

        self.assertEqual(counts['failed'], 1)
        self.assertEqual(counts['completed'], 1)

        work_pmid = Path(workdir) / '99999'
        self.assertTrue((work_pmid / 'GCST_BAD.failed.yaml').exists())
        self.assertTrue((work_pmid / 'GCST_GOOD.tsv.gz').exists())

    def test_failed_yaml_contains_error_and_timestamp(self):
        """failed.yaml has gcst_id, error, timestamp."""
        path = Path(self.tmpdir) / 'test.failed.yaml'
        _write_failed_yaml(path, 'GCST_X', 'Something went wrong')
        with open(path) as fh:
            data = yaml.safe_load(fh)
        self.assertEqual(data['gcst_id'], 'GCST_X')
        self.assertIn('error', data)
        self.assertIn('timestamp', data)
        self.assertIn('Something went wrong', data['error'])

    def test_summary_counts(self):
        """build() returns dict with completed/skipped/failed counts."""
        disc, workdir = self._make_env()
        counts = Stage1Builder().build(disc, workdir=workdir)
        self.assertIn('completed', counts)
        self.assertIn('skipped', counts)
        self.assertIn('failed', counts)


if __name__ == '__main__':
    unittest.main(verbosity=2)
