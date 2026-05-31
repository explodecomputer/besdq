"""Tests for the distance-based peak counter."""

import unittest
from besdq.gwas_ssf_reader import GwasSsfRow
from besdq.peak_count import count_peaks_by_distance


def _row(chr: str, bp: int) -> GwasSsfRow:
    return GwasSsfRow(chr=chr, bp=bp, a1='A', a2='T',
                      rsid=None, beta=0.1, se=0.01, eaf=0.3, p=1e-5)


class TestCountPeaksByDistance(unittest.TestCase):

    def test_empty_returns_zero(self):
        self.assertEqual(count_peaks_by_distance([], radius=500_000), 0)

    def test_single_row_is_one_peak(self):
        self.assertEqual(count_peaks_by_distance([_row('1', 1_000_000)], 500_000), 1)

    def test_all_within_radius_is_one_peak(self):
        rows = [_row('1', 1_000_000), _row('1', 1_200_000), _row('1', 1_400_000)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 1)

    def test_all_beyond_radius_each_independent(self):
        rows = [_row('1', 1_000_000), _row('1', 2_000_000), _row('1', 3_000_000)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 3)

    def test_exact_boundary_within_is_same_peak(self):
        # distance == radius → still same peak (not strictly greater)
        rows = [_row('1', 1_000_000), _row('1', 1_500_000)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 1)

    def test_one_beyond_boundary_is_new_peak(self):
        rows = [_row('1', 1_000_000), _row('1', 1_500_001)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 2)

    def test_different_chromosomes_each_independent(self):
        rows = [_row('1', 1_000_000), _row('2', 1_000_000), _row('3', 1_000_000)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 3)

    def test_mixed_chromosomes_grouped_correctly(self):
        rows = [
            _row('1', 1_000_000), _row('1', 1_200_000),  # same peak on chr1
            _row('2', 1_000_000), _row('2', 2_000_000),  # two peaks on chr2
        ]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 3)

    def test_unsorted_input_handled(self):
        rows = [_row('1', 3_000_000), _row('1', 1_000_000), _row('1', 2_000_000)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 3)

    def test_alpha_chromosome_names(self):
        rows = [_row('X', 1_000_000), _row('X', 1_200_000), _row('Y', 1_000_000)]
        self.assertEqual(count_peaks_by_distance(rows, radius=500_000), 2)
