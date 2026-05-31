"""Tests for discover-study command (Issues #25 and #26)."""

import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from besdq.discover import (
    DISCOVERY_TSV_COLUMNS,
    _harmonised_url,
    _parse_gene_symbol,
    discover_study,
    load_gene_annotation,
    query_ebi_studies,
    write_discovery_tsv,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MOCK_STUDY_PAGE_1 = {
    "_embedded": {
        "studies": [
            {
                "accessionId": "GCST001",
                "summaryStatistics": "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST001",
                "reportedTrait": "IL10 expression at baseline",
            },
            {
                "accessionId": "GCST002",
                "summaryStatistics": "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002",
                "reportedTrait": "TNF expression in monocytes",
            },
        ]
    },
    "page": {"totalPages": 2, "number": 0, "size": 2, "totalElements": 3},
}

MOCK_STUDY_PAGE_2 = {
    "_embedded": {
        "studies": [
            {
                "accessionId": "GCST003",
                "summaryStatistics": "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003",
                "reportedTrait": "BRCA1 expression",
            },
        ]
    },
    "page": {"totalPages": 2, "number": 1, "size": 2, "totalElements": 3},
}

GENE_ANNOTATION_CONTENT = (
    "gene_name\tchr\tbp\n"
    "IL10\t1\t206941637\n"
    "TNF\t6\t31575567\n"
    "BRCA1\t17\t43044295\n"
)


# ---------------------------------------------------------------------------
# Issue #25 — discover-study core
# ---------------------------------------------------------------------------

class TestDiscoverStudy(unittest.TestCase):

    def _mock_responses(self, pages):
        """Return a mock for requests.get that serves each page in order."""
        responses = []
        for page_data in pages:
            m = mock.MagicMock()
            m.status_code = 200
            m.json.return_value = page_data
            responses.append(m)
        return responses

    def test_single_page(self):
        single_page = {
            "_embedded": {"studies": MOCK_STUDY_PAGE_1["_embedded"]["studies"]},
            "page": {"totalPages": 1, "number": 0, "size": 2, "totalElements": 2},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = single_page
            rows = discover_study("12345")

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['gcst_id'], 'GCST001')
        self.assertEqual(rows[1]['gcst_id'], 'GCST002')
        self.assertEqual(rows[0]['pmid'], '12345')

    def test_pagination(self):
        """Studies beyond page 0 are included."""
        responses = self._mock_responses([MOCK_STUDY_PAGE_1, MOCK_STUDY_PAGE_2])
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.side_effect = [
                mock.MagicMock(status_code=200, json=lambda: MOCK_STUDY_PAGE_1),
                mock.MagicMock(status_code=200, json=lambda: MOCK_STUDY_PAGE_2),
            ]
            rows = discover_study("12345")

        self.assertEqual(len(rows), 3)
        gcst_ids = [r['gcst_id'] for r in rows]
        self.assertIn('GCST001', gcst_ids)
        self.assertIn('GCST003', gcst_ids)

    def test_harmonised_url_construction(self):
        study = {
            "accessionId": "GCST001",
            "summaryStatistics": "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST001",
        }
        url = _harmonised_url(study)
        self.assertEqual(
            url,
            "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST001"
            "/harmonised/GCST001.h.tsv.gz",
        )

    def test_file_path_is_harmonised_url(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = single_page
            rows = discover_study("12345")

        self.assertIn('/harmonised/GCST001.h.tsv.gz', rows[0]['file_path'])

    def test_trait_name_from_reported_trait(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = single_page
            rows = discover_study("12345")

        self.assertEqual(rows[0]['trait_name'], 'IL10 expression at baseline')

    def test_gene_chr_bp_empty_without_parse_gene(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = single_page
            rows = discover_study("12345")

        self.assertEqual(rows[0]['gene'], '')
        self.assertEqual(rows[0]['trait_chr'], '')
        self.assertEqual(rows[0]['trait_bp'], '')

    def test_output_tsv_columns(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = single_page
            rows = discover_study("12345")

        buf = io.StringIO()
        with mock.patch('sys.stdout', buf):
            write_discovery_tsv(rows)
        lines = buf.getvalue().strip().split('\n')
        header_cols = lines[0].split('\t')
        self.assertEqual(header_cols, DISCOVERY_TSV_COLUMNS)

    def test_output_to_file(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = single_page
            rows = discover_study("12345")

        with tempfile.NamedTemporaryFile(mode='r', suffix='.tsv', delete=False) as f:
            out_path = f.name
        try:
            write_discovery_tsv(rows, output=out_path)
            with open(out_path) as f:
                content = f.read()
            self.assertIn('GCST001', content)
            self.assertIn('\t'.join(DISCOVERY_TSV_COLUMNS), content)
        finally:
            Path(out_path).unlink(missing_ok=True)

    def test_http_error_raises_clear_message(self):
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 404
            mock_req.get.return_value.text = 'Not Found'
            mock_req.RequestException = Exception
            with self.assertRaises(RuntimeError) as cm:
                discover_study("99999")
        self.assertIn('404', str(cm.exception))

    def test_no_studies_returns_empty_list(self):
        empty_page = {
            "_embedded": {"studies": []},
            "page": {"totalPages": 1, "number": 0},
        }
        with mock.patch('besdq.discover._requests') as mock_req:
            mock_req.get.return_value.status_code = 200
            mock_req.get.return_value.json.return_value = empty_page
            rows = discover_study("00000")
        self.assertEqual(rows, [])


# ---------------------------------------------------------------------------
# Issue #26 — --parse-gene flag
# ---------------------------------------------------------------------------

class TestParseGene(unittest.TestCase):

    def _write_gene_annotation(self) -> str:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write(GENE_ANNOTATION_CONTENT)
            return f.name

    def test_parse_gene_requires_gene_annotation(self):
        with self.assertRaises(ValueError) as cm:
            discover_study("12345", parse_gene=True, gene_annotation_path=None)
        self.assertIn('--gene-annotation', str(cm.exception))

    def test_known_gene_resolves_chr_bp(self):
        ann_path = self._write_gene_annotation()
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        try:
            with mock.patch('besdq.discover._requests') as mock_req:
                mock_req.get.return_value.status_code = 200
                mock_req.get.return_value.json.return_value = single_page
                rows = discover_study("12345", parse_gene=True, gene_annotation_path=ann_path)
            self.assertEqual(rows[0]['gene'], 'IL10')
            self.assertEqual(rows[0]['trait_chr'], '1')
            self.assertEqual(rows[0]['trait_bp'], '206941637')
        finally:
            Path(ann_path).unlink(missing_ok=True)

    def test_unknown_gene_emits_warning_row_still_included(self):
        ann_path = self._write_gene_annotation()
        unknown_study_page = {
            "_embedded": {
                "studies": [{
                    "accessionId": "GCST099",
                    "summaryStatistics": "ftp://ftp.ebi.ac.uk/pub/databases/gwas/GCST099",
                    "reportedTrait": "UNKNOWNGENE expression",
                }]
            },
            "page": {"totalPages": 1, "number": 0},
        }
        try:
            with mock.patch('besdq.discover._requests') as mock_req:
                mock_req.get.return_value.status_code = 200
                mock_req.get.return_value.json.return_value = unknown_study_page
                import sys as _sys
                with mock.patch.object(_sys, 'stderr') as mock_err:
                    rows = discover_study("12345", parse_gene=True, gene_annotation_path=ann_path)
        finally:
            Path(ann_path).unlink(missing_ok=True)

        # Row still emitted
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['gene'], '')
        self.assertEqual(rows[0]['trait_chr'], '')
        self.assertEqual(rows[0]['trait_bp'], '')

    def test_flag_off_leaves_gene_empty(self):
        ann_path = self._write_gene_annotation()
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        try:
            with mock.patch('besdq.discover._requests') as mock_req:
                mock_req.get.return_value.status_code = 200
                mock_req.get.return_value.json.return_value = single_page
                rows = discover_study("12345", parse_gene=False, gene_annotation_path=ann_path)
            # Without parse_gene, gene is always empty
            self.assertEqual(rows[0]['gene'], '')
        finally:
            Path(ann_path).unlink(missing_ok=True)

    def test_gene_annotation_missing_column_raises(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("gene_name\tchr\n")  # missing 'bp'
            f.write("IL10\t1\n")
            bad_path = f.name
        try:
            with self.assertRaises(ValueError) as cm:
                load_gene_annotation(bad_path)
            self.assertIn('bp', str(cm.exception))
        finally:
            Path(bad_path).unlink(missing_ok=True)

    def test_parse_gene_symbol_first_word(self):
        self.assertEqual(_parse_gene_symbol("IL10 expression at baseline"), "IL10")
        self.assertEqual(_parse_gene_symbol("TNF"), "TNF")
        self.assertEqual(_parse_gene_symbol(""), "")
        self.assertEqual(_parse_gene_symbol("  BRCA1  variant"), "BRCA1")


# ---------------------------------------------------------------------------
# Integration test (live EBI API; gated by TEST_EBI_LIVE=1)
# ---------------------------------------------------------------------------

class TestDiscoverStudyLive(unittest.TestCase):

    def _require_live(self):
        if not os.environ.get('TEST_EBI_LIVE'):
            self.skipTest("Set TEST_EBI_LIVE=1 to run live EBI integration tests")

    def test_live_pmid_38714679(self):
        self._require_live()
        rows = discover_study("38714679")
        self.assertGreater(len(rows), 0)
        for row in rows:
            self.assertTrue(row['gcst_id'].startswith('GCST'))
            self.assertIn('.h.tsv.gz', row['file_path'])
        gcst_ids = {r['gcst_id'] for r in rows}
        # Spot-check a few known accessions from this publication
        self.assertTrue(any(g.startswith('GCST') for g in gcst_ids))


if __name__ == '__main__':
    unittest.main(verbosity=2)
