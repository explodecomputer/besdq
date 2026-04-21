"""Unit tests for BESD query functionality."""

import unittest
import tempfile
from pathlib import Path
from besdq import BESDQueryEngine, BESDIndexBuilder, BESDQueryIndex


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


if __name__ == '__main__':
    unittest.main(verbosity=2)



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


if __name__ == '__main__':
    unittest.main(verbosity=2)
