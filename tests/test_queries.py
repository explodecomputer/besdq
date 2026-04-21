"""Unit tests for BESD query functionality."""

import unittest
import tempfile
from pathlib import Path
from unittest import mock
import io
import contextlib
from besdq import BESDQueryEngine, BESDIndexBuilder, BESDQueryIndex
from besdq.besd_reader import IndexReader, calculate_p_value as reader_calculate_p_value
from besdq.sqlite_query import calculate_p_value as sqlite_calculate_p_value
from besdq import cli



# Path to test data
DATA_DIR = Path(__file__).parent.parent / "data"
WESTRA_BESD = str(DATA_DIR / "westra_eqtl_hg19")


class TestBESDQueryEngine(unittest.TestCase):
    """Test BESD query engine with westra_eqtl_hg19 dataset."""

    @classmethod
    def setUpClass(cls):
        """Initialize query engine with westra data."""
        cls.engine = BESDQueryEngine(WESTRA_BESD)

    def test_load_data(self):
        """Test that BESD data loads correctly."""
        self.assertEqual(len(self.engine.snps), 506049)
        self.assertEqual(len(self.engine.probes), 5966)
        self.assertEqual(self.engine.besd.format_type, '3F')

    def test_single_position_query(self):
        """Test querying at exact SNP and probe positions."""
        # Query at exact SNP and probe positions from earlier testing
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        self.assertEqual(associations[0]['snp_id'], 'rs3818646')
        self.assertEqual(associations[0]['probe_id'], 'ILMN_2349633')

    def test_range_query(self):
        """Test querying a range of SNP and probe positions."""
        # Query from earlier testing: 1:100kb-2000kb SNPs, 1:1000kb-2000kb probes
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        self.assertEqual(len(associations), 175)

    def test_p_value_calculation(self):
        """Test that p-values are calculated correctly."""
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        assoc = associations[0]

        # From earlier testing, we know:
        # beta = -0.436080, se = 0.040022
        self.assertAlmostEqual(assoc['beta'], -0.436080, places=4)
        self.assertAlmostEqual(assoc['se'], 0.040022, places=4)
        self.assertGreaterEqual(assoc['pval'], 0)  # p-value should be non-negative

    def test_beta_and_se_stored(self):
        """Test that betas and SEs are stored correctly."""
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        # Verify all associations have beta and SE values
        for assoc in associations:
            self.assertIn('beta', assoc)
            self.assertIn('se', assoc)
            self.assertIsInstance(assoc['beta'], float)
            self.assertIsInstance(assoc['se'], float)
            self.assertGreater(assoc['se'], 0)  # SEs should be positive

    def test_snp_indexing(self):
        """Test that SNP metadata is correctly indexed."""
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        assoc = associations[0]

        # Verify SNP fields
        self.assertEqual(assoc['snp_id'], 'rs3818646')
        self.assertEqual(assoc['snp_chr'], '1')
        self.assertEqual(assoc['snp_bp'], 1191870)
        self.assertEqual(assoc['a1'], 'T')
        self.assertEqual(assoc['a2'], 'C')

    def test_probe_indexing(self):
        """Test that probe metadata is correctly indexed."""
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        assoc = associations[0]

        # Verify probe fields
        self.assertEqual(assoc['probe_id'], 'ILMN_2349633')
        self.assertEqual(assoc['probe_chr'], '1')
        self.assertEqual(assoc['probe_bp'], 1140818)
        self.assertEqual(assoc['gene'], 'TNFRSF18')

    def test_chromosome_filtering(self):
        """Test that queries respect chromosome boundaries."""
        # Query for chr 1
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )

        # All results should be from chr 1
        for assoc in associations:
            self.assertEqual(assoc['snp_chr'], '1')
            self.assertEqual(assoc['probe_chr'], '1')

    def test_position_filtering(self):
        """Test that queries respect position boundaries."""
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )

        # All results should be within the queried ranges
        for assoc in associations:
            self.assertGreaterEqual(assoc['snp_bp'], 100000)
            self.assertLessEqual(assoc['snp_bp'], 2000000)
            self.assertGreaterEqual(assoc['probe_bp'], 1000000)
            self.assertLessEqual(assoc['probe_bp'], 2000000)

    def test_empty_query(self):
        """Test that queries with no matches return empty list."""
        # Query a range with no associations
        associations = self.engine.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=110,  # Small range
            probe_chr='1', probe_start_kb=1900, probe_end_kb=2000,  # Far from SNPs
        )
        self.assertEqual(len(associations), 0)

    def test_multi_snp_query(self):
        """Test querying multiple SNPs."""
        # Test query_by_snp_id with multiple SNPs
        assocs1 = self.engine.query_by_snp_id('rs3818646')
        assocs2 = self.engine.query_by_snp_id('rs7515488')

        self.assertEqual(len(assocs1), 5)
        self.assertEqual(len(assocs2), 6)

        # Combined
        combined = assocs1 + assocs2
        self.assertEqual(len(combined), 11)

    def test_multi_probe_query(self):
        """Test querying multiple probes."""
        # Test query_by_probe_id with multiple probes
        assocs1 = self.engine.query_by_probe_id('ILMN_2349633')
        assocs2 = self.engine.query_by_probe_id('ILMN_2112256')

        self.assertEqual(len(assocs1), 20)
        self.assertEqual(len(assocs2), 20)

        # Check they're different probes
        probe_ids = {a['probe_id'] for a in assocs1 + assocs2}
        self.assertEqual(len(probe_ids), 2)

    def test_gene_query(self):
        """Test querying by gene name."""
        associations = self.engine.query_by_gene('TNFRSF18')
        self.assertEqual(len(associations), 20)

        # All results should be for the same gene
        for assoc in associations:
            self.assertEqual(assoc['gene'], 'TNFRSF18')


class TestBESDQueryIndex(unittest.TestCase):
    """Test SQLite-indexed BESD query engine."""

    @classmethod
    def setUpClass(cls):
        """Create SQLite index database for testing."""
        cls.temp_db = tempfile.NamedTemporaryFile(suffix='.db', delete=False)
        cls.db_path = cls.temp_db.name
        cls.temp_db.close()

        # Build index
        builder = BESDIndexBuilder(cls.db_path)
        builder.build(WESTRA_BESD, force=True)

        # Initialize query engine
        cls.index = BESDQueryIndex(cls.db_path)

    @classmethod
    def tearDownClass(cls):
        """Clean up temporary database."""
        cls.index.close()
        Path(cls.db_path).unlink()

    def test_metadata_loading(self):
        """Test that metadata loads correctly."""
        self.assertEqual(self.index.metadata['format_type'], '3F')
        self.assertEqual(self.index.metadata['n_snps'], '506049')
        self.assertEqual(self.index.metadata['n_probes'], '5966')

    def test_single_position_query(self):
        """Test querying at exact SNP and probe positions."""
        associations = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        self.assertEqual(associations[0]['snp_id'], 'rs3818646')
        self.assertEqual(associations[0]['probe_id'], 'ILMN_2349633')

    def test_range_query(self):
        """Test querying a range of SNP and probe positions."""
        associations = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=100, snp_end_kb=2000,
            probe_chr='1', probe_start_kb=1000, probe_end_kb=2000,
        )
        self.assertEqual(len(associations), 175)

    def test_p_value_calculation(self):
        """Test that p-values are calculated correctly."""
        associations = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        self.assertEqual(len(associations), 1)
        assoc = associations[0]

        self.assertAlmostEqual(assoc['beta'], -0.436080, places=4)
        self.assertAlmostEqual(assoc['se'], 0.040022, places=4)
        self.assertGreaterEqual(assoc['pval'], 0)

    def test_query_by_probe_id(self):
        """Test querying by probe ID."""
        associations = self.index.query_by_probe_id('ILMN_2349633')
        self.assertEqual(len(associations), 20)
        self.assertTrue(all(a['probe_id'] == 'ILMN_2349633' for a in associations))

    def test_query_by_snp_id(self):
        """Test querying by SNP ID."""
        associations = self.index.query_by_snp_id('rs3818646')
        self.assertEqual(len(associations), 5)
        self.assertTrue(all(a['snp_id'] == 'rs3818646' for a in associations))

    def test_consistency_with_besd_reader(self):
        """Test that SQLite index produces same results as BESD reader."""
        engine = BESDQueryEngine(WESTRA_BESD)

        # Compare results for single position query
        besd_results = engine.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )
        index_results = self.index.query_cis_window(
            snp_chr='1', snp_start_kb=1191.87, snp_end_kb=1191.87,
            probe_chr='1', probe_start_kb=1140.818, probe_end_kb=1140.818,
        )

        self.assertEqual(len(besd_results), len(index_results))
        self.assertEqual(besd_results[0]['snp_id'], index_results[0]['snp_id'])
        self.assertEqual(besd_results[0]['probe_id'], index_results[0]['probe_id'])
        self.assertAlmostEqual(besd_results[0]['beta'], index_results[0]['beta'], places=5)
        self.assertAlmostEqual(besd_results[0]['se'], index_results[0]['se'], places=5)

    def test_multi_snp_query_index(self):
        """Test querying multiple SNPs from index."""
        assocs1 = self.index.query_by_snp_id('rs3818646')
        assocs2 = self.index.query_by_snp_id('rs7515488')

        self.assertEqual(len(assocs1), 5)
        self.assertEqual(len(assocs2), 6)

    def test_multi_probe_query_index(self):
        """Test querying multiple probes from index."""
        assocs1 = self.index.query_by_probe_id('ILMN_2349633')
        assocs2 = self.index.query_by_probe_id('ILMN_2112256')

        self.assertEqual(len(assocs1), 20)
        self.assertEqual(len(assocs2), 20)

        # Check they're different probes
        probe_ids = {a['probe_id'] for a in assocs1 + assocs2}
        self.assertEqual(len(probe_ids), 2)

    def test_gene_query_index(self):
        """Test querying by gene name from index."""
        associations = self.index.query_by_gene('TNFRSF18')
        self.assertEqual(len(associations), 20)

        # All results should be for the same gene
        for assoc in associations:
            self.assertEqual(assoc['gene'], 'TNFRSF18')

    def test_consistency_snp_query(self):
        """Test that SNP queries are consistent between BESD and index."""
        engine = BESDQueryEngine(WESTRA_BESD)

        besd_results = engine.query_by_snp_id('rs3818646')
        index_results = self.index.query_by_snp_id('rs3818646')

        self.assertEqual(len(besd_results), len(index_results))
        # Check same associations returned
        besd_assoc_ids = {(a['snp_id'], a['probe_id']) for a in besd_results}
        index_assoc_ids = {(a['snp_id'], a['probe_id']) for a in index_results}
        self.assertEqual(besd_assoc_ids, index_assoc_ids)

    def test_consistency_probe_query(self):
        """Test that probe queries are consistent between BESD and index."""
        engine = BESDQueryEngine(WESTRA_BESD)

        besd_results = engine.query_by_probe_id('ILMN_2349633')
        index_results = self.index.query_by_probe_id('ILMN_2349633')

        self.assertEqual(len(besd_results), len(index_results))
        # Check same associations returned
        besd_assoc_ids = {(a['snp_id'], a['probe_id']) for a in besd_results}
        index_assoc_ids = {(a['snp_id'], a['probe_id']) for a in index_results}
        self.assertEqual(besd_assoc_ids, index_assoc_ids)

    def test_consistency_gene_query(self):
        """Test that gene queries are consistent between BESD and index."""
        engine = BESDQueryEngine(WESTRA_BESD)

        besd_results = engine.query_by_gene('TNFRSF18')
        index_results = self.index.query_by_gene('TNFRSF18')

        self.assertEqual(len(besd_results), len(index_results))
        # Check same associations returned
        besd_assoc_ids = {(a['snp_id'], a['probe_id']) for a in besd_results}
        index_assoc_ids = {(a['snp_id'], a['probe_id']) for a in index_results}
        self.assertEqual(besd_assoc_ids, index_assoc_ids)


class TestIndexReaderParsing(unittest.TestCase):
    """Test index file parsing edge cases."""

    def test_read_esi_row_idx_ignores_non_data_lines(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.esi', delete=False) as f:
            f.write("# comment\n")
            f.write("\n")
            f.write("1 rsA 0.1 100 A G 0.4\n")
            f.write("bad line\n")
            f.write("1 rsB 0.2 200 C T 0.2\n")
            esi_path = f.name

        try:
            snps = IndexReader.read_esi(esi_path)
            self.assertEqual(len(snps), 2)
            self.assertEqual(snps[0]['row_idx'], 0)
            self.assertEqual(snps[1]['row_idx'], 1)
            self.assertEqual(snps[0]['snp_id'], 'rsA')
            self.assertEqual(snps[1]['snp_id'], 'rsB')
        finally:
            Path(esi_path).unlink()

    def test_read_epi_row_idx_ignores_non_data_lines(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.epi', delete=False) as f:
            f.write("# comment\n")
            f.write("1 probeA 0.0 1000 GENE1 +\n")
            f.write("\n")
            f.write("bad line\n")
            f.write("1 probeB 0.0 2000 GENE2 -\n")
            epi_path = f.name

        try:
            probes = IndexReader.read_epi(epi_path)
            self.assertEqual(len(probes), 2)
            self.assertEqual(probes[0]['row_idx'], 0)
            self.assertEqual(probes[1]['row_idx'], 1)
            self.assertEqual(probes[0]['probe_id'], 'probeA')
            self.assertEqual(probes[1]['probe_id'], 'probeB')
        finally:
            Path(epi_path).unlink()


class TestPValueEdgeCases(unittest.TestCase):
    """Test numeric stability of p-value computations."""

    def test_p_value_extreme_z_is_bounded_reader(self):
        pval = reader_calculate_p_value(beta=1e9, se=1e-12)
        self.assertGreaterEqual(pval, 0.0)
        self.assertLessEqual(pval, 1.0)

    def test_p_value_extreme_z_is_bounded_sqlite(self):
        pval = sqlite_calculate_p_value(beta=1e9, se=1e-12)
        self.assertGreaterEqual(pval, 0.0)
        self.assertLessEqual(pval, 1.0)

    def test_p_value_zero_or_negative_se_defaults_to_one(self):
        self.assertEqual(reader_calculate_p_value(beta=1.0, se=0.0), 1.0)
        self.assertEqual(reader_calculate_p_value(beta=1.0, se=-1.0), 1.0)
        self.assertEqual(sqlite_calculate_p_value(beta=1.0, se=0.0), 1.0)
        self.assertEqual(sqlite_calculate_p_value(beta=1.0, se=-1.0), 1.0)


class TestCLIValidation(unittest.TestCase):
    """Test CLI parsing and argument validation edge cases."""

    def test_parse_chrpos_reversed_range_raises(self):
        with self.assertRaises(ValueError):
            cli.parse_chrpos("1:200-100")

    def test_parse_chrpos_malformed_position_raises(self):
        with self.assertRaises(ValueError):
            cli.parse_chrpos("1:abc")

    def test_cli_rejects_conflicting_identifier_queries(self):
        test_argv = [
            "besdq",
            "--beqtl-summary", "dummy",
            "--out", "out",
            "--snp", "rs1",
            "--probe", "probe1",
        ]
        with (
            mock.patch("sys.argv", test_argv),
            mock.patch("besdq.cli.BESDQueryEngine"),
            mock.patch("besdq.cli.BESDQueryIndex"),
            contextlib.redirect_stderr(io.StringIO()) as stderr
        ):
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("mutually exclusive", stderr.getvalue())

    def test_cli_rejects_mixed_identifier_and_region_queries(self):
        test_argv = [
            "besdq",
            "--beqtl-summary", "dummy",
            "--out", "out",
            "--snp", "rs1",
            "--snp-chrpos", "1:100-200",
            "--probe-chrpos", "1:100-200",
        ]
        with (
            mock.patch("sys.argv", test_argv),
            mock.patch("besdq.cli.BESDQueryEngine"),
            mock.patch("besdq.cli.BESDQueryIndex"),
            contextlib.redirect_stderr(io.StringIO()) as stderr
        ):
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("Cannot combine --snp/--probe/--gene", stderr.getvalue())


if __name__ == '__main__':
    unittest.main(verbosity=2)
