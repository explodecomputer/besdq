"""Tests for discover-study command (Issues #25 and #26)."""

import io
import os
import sys
import tempfile
import unittest
from contextlib import contextmanager
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
    "gene_name\tchr\tbp\tcanonical\n"
    "IL10\t1\t206941637\t\n"
    "TNF\t6\t31575567\t\n"
    "BRCA1\t17\t43044295\t\n"
    "IL1RA\t2\t113117981\tIL1RN\n"
    "TNFA\t6\t31575567\tTNF\n"
    "IL1RN\t2\t113117981\t\n"
)


# ---------------------------------------------------------------------------
# Mock helpers
# ---------------------------------------------------------------------------

def _api_response(page_data: dict):
    m = mock.MagicMock()
    m.status_code = 200
    m.json.return_value = page_data
    return m


def _meta_yaml_response(sample_size: int = 1060):
    m = mock.MagicMock()
    m.status_code = 200
    m.text = f"samples:\n- sample_size: {sample_size}\n"
    return m


def _meta_yaml_404():
    return mock.MagicMock(status_code=404, text="Not Found")


def _side_effects_for_page(page_data: dict, sample_size: int = 1060):
    """API page response + one meta YAML per study.

    Valid for single-page calls where all API pages are fetched before
    _fetch_sample_size is called. For multi-page tests, build the sequence
    manually (all API pages first, then all meta YAMLs).
    """
    n = len(page_data["_embedded"]["studies"])
    return [_api_response(page_data)] + [_meta_yaml_response(sample_size)] * n


@contextmanager
def patch_requests(side_effects):
    """Patch _requests and _REQUESTS_AVAILABLE in besdq.discover together."""
    with mock.patch('besdq.discover._requests') as mock_req, \
         mock.patch('besdq.discover._REQUESTS_AVAILABLE', True):
        mock_req.RequestException = Exception
        mock_req.get.side_effect = list(side_effects)
        yield mock_req


# ---------------------------------------------------------------------------
# Issue #25 — discover-study core
# ---------------------------------------------------------------------------

class TestDiscoverStudy(unittest.TestCase):

    def test_single_page(self):
        single_page = {
            "_embedded": {"studies": MOCK_STUDY_PAGE_1["_embedded"]["studies"]},
            "page": {"totalPages": 1, "number": 0, "size": 2, "totalElements": 2},
        }
        with patch_requests(_side_effects_for_page(single_page)):
            rows = discover_study("12345")

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['gcst_id'], 'GCST001')
        self.assertEqual(rows[1]['gcst_id'], 'GCST002')
        self.assertEqual(rows[0]['pmid'], '12345')

    def test_pagination(self):
        """Studies beyond page 0 are included."""
        # query_ebi_studies fetches all API pages first, then _fetch_sample_size
        # is called per study — so API responses come before meta YAML responses.
        all_studies = (
            MOCK_STUDY_PAGE_1["_embedded"]["studies"] +
            MOCK_STUDY_PAGE_2["_embedded"]["studies"]
        )
        side_effects = (
            [_api_response(MOCK_STUDY_PAGE_1), _api_response(MOCK_STUDY_PAGE_2)] +
            [_meta_yaml_response() for _ in all_studies]
        )
        with patch_requests(side_effects):
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
        with patch_requests(_side_effects_for_page(single_page)):
            rows = discover_study("12345")
        self.assertIn('/harmonised/GCST001.h.tsv.gz', rows[0]['file_path'])

    def test_trait_name_from_reported_trait(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests(_side_effects_for_page(single_page)):
            rows = discover_study("12345")
        self.assertEqual(rows[0]['trait_name'], 'IL10 expression at baseline')

    def test_gene_chr_bp_empty_without_parse_gene(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests(_side_effects_for_page(single_page)):
            rows = discover_study("12345")
        self.assertEqual(rows[0]['gene'], '')
        self.assertEqual(rows[0]['trait_chr'], '')
        self.assertEqual(rows[0]['trait_bp'], '')

    def test_output_tsv_columns(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests(_side_effects_for_page(single_page)):
            rows = discover_study("12345")

        buf = io.StringIO()
        with mock.patch('sys.stdout', buf):
            write_discovery_tsv(rows)
        lines = buf.getvalue().strip().split('\n')
        self.assertEqual(lines[0].split('\t'), DISCOVERY_TSV_COLUMNS)

    def test_output_to_file(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests(_side_effects_for_page(single_page)):
            rows = discover_study("12345")

        with tempfile.NamedTemporaryFile(mode='r', suffix='.tsv', delete=False) as f:
            out_path = f.name
        try:
            write_discovery_tsv(rows, output=out_path)
            content = Path(out_path).read_text()
            self.assertIn('GCST001', content)
            self.assertIn('\t'.join(DISCOVERY_TSV_COLUMNS), content)
        finally:
            Path(out_path).unlink(missing_ok=True)

    def test_http_error_raises_clear_message(self):
        error_resp = mock.MagicMock(status_code=404, text='Not Found')
        with patch_requests([error_resp]):
            with self.assertRaises(RuntimeError) as cm:
                discover_study("99999")
        self.assertIn('404', str(cm.exception))

    def test_sample_size_populated_from_meta_yaml(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests(_side_effects_for_page(single_page, sample_size=575)):
            rows = discover_study("12345")
        self.assertEqual(rows[0]['sample_size'], 575)

    def test_sample_size_empty_when_meta_yaml_unavailable(self):
        single_page = {
            "_embedded": {"studies": [MOCK_STUDY_PAGE_1["_embedded"]["studies"][0]]},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests([_api_response(single_page), _meta_yaml_404()]):
            rows = discover_study("12345")
        self.assertEqual(rows[0]['sample_size'], '')

    def test_no_studies_returns_empty_list(self):
        empty_page = {
            "_embedded": {"studies": []},
            "page": {"totalPages": 1, "number": 0},
        }
        with patch_requests([_api_response(empty_page)]):
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
            with patch_requests(_side_effects_for_page(single_page)):
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
            with patch_requests(_side_effects_for_page(unknown_study_page)):
                rows = discover_study("12345", parse_gene=True, gene_annotation_path=ann_path)
        finally:
            Path(ann_path).unlink(missing_ok=True)

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
            with patch_requests(_side_effects_for_page(single_page)):
                rows = discover_study("12345", parse_gene=False, gene_annotation_path=ann_path)
            self.assertEqual(rows[0]['gene'], '')
        finally:
            Path(ann_path).unlink(missing_ok=True)

    def test_gene_annotation_missing_column_raises(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("gene_name\tchr\n")
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

    def test_alias_resolves_to_canonical_symbol(self):
        ann_path = self._write_gene_annotation()
        il1ra_page = {
            "_embedded": {
                "studies": [{
                    "accessionId": "GCST_IL1",
                    "summaryStatistics": "ftp://ftp.ebi.ac.uk/pub/databases/gwas/GCST_IL1",
                    "reportedTrait": "IL1Ra expression at baseline",
                }]
            },
            "page": {"totalPages": 1, "number": 0},
        }
        try:
            with patch_requests(_side_effects_for_page(il1ra_page)):
                rows = discover_study("12345", parse_gene=True, gene_annotation_path=ann_path)
            self.assertEqual(rows[0]['gene'], 'IL1RN')
            self.assertEqual(rows[0]['trait_chr'], '2')
            self.assertEqual(rows[0]['trait_bp'], '113117981')
        finally:
            Path(ann_path).unlink(missing_ok=True)

    def test_case_insensitive_lookup(self):
        ann_path = self._write_gene_annotation()
        try:
            for variant in ('IL10', 'il10', 'Il10'):
                page = {
                    "_embedded": {
                        "studies": [{
                            "accessionId": "GCST_CI",
                            "summaryStatistics": "ftp://ftp.ebi.ac.uk/s",
                            "reportedTrait": f"{variant} expression",
                        }]
                    },
                    "page": {"totalPages": 1, "number": 0},
                }
                with patch_requests(_side_effects_for_page(page)):
                    rows = discover_study("12345", parse_gene=True, gene_annotation_path=ann_path)
                self.assertEqual(rows[0]['trait_chr'], '1', f"Failed for variant {variant!r}")
        finally:
            Path(ann_path).unlink(missing_ok=True)

    def test_previous_symbol_resolves_to_canonical(self):
        ann_path = self._write_gene_annotation()
        page = {
            "_embedded": {
                "studies": [{
                    "accessionId": "GCST_TNF",
                    "summaryStatistics": "ftp://ftp.ebi.ac.uk/s",
                    "reportedTrait": "TNFA expression",
                }]
            },
            "page": {"totalPages": 1, "number": 0},
        }
        try:
            with patch_requests(_side_effects_for_page(page)):
                rows = discover_study("12345", parse_gene=True, gene_annotation_path=ann_path)
            self.assertEqual(rows[0]['gene'], 'TNF')
            self.assertEqual(rows[0]['trait_chr'], '6')
        finally:
            Path(ann_path).unlink(missing_ok=True)


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
            self.assertIsInstance(row['sample_size'], int)
            self.assertGreater(row['sample_size'], 0)
        gcst_ids = {r['gcst_id'] for r in rows}
        self.assertTrue(any(g.startswith('GCST') for g in gcst_ids))


if __name__ == '__main__':
    unittest.main(verbosity=2)
