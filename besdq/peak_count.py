"""Distance-based independent peak counter.

Used as a fallback when plink2 LD clumping was not run: counts variants that
are more than `radius` bp from the nearest already-counted peak on the same
chromosome as independent peaks (greedy positional sweep).

This approximation can over-count in high-LD regions.  Callers should flag
the result as an estimate when this function is used instead of clumping.
"""

from typing import List

from .gwas_ssf_reader import GwasSsfRow


def count_peaks_by_distance(rows: List[GwasSsfRow], radius: int) -> int:
    """Count independent peaks by greedy positional sweep.

    Sorts `rows` by (chromosome, bp) then greedily assigns each row to a new
    peak if it is more than `radius` bp from the most recent peak on the same
    chromosome.

    Args:
        rows:   SNP rows to count peaks in (any order).
        radius: Minimum bp distance between independent peaks.

    Returns:
        Number of independent peaks (0 if rows is empty).
    """
    if not rows:
        return 0

    # Sort by chromosome (numeric then alpha) then position
    def _sort_key(r: GwasSsfRow):
        try:
            return (0, int(r.chr), r.bp)
        except ValueError:
            return (1, r.chr, r.bp)

    sorted_rows = sorted(rows, key=_sort_key)

    peaks = 0
    last_chr: str = ""
    last_bp: int = -1

    for row in sorted_rows:
        if row.chr != last_chr or abs(row.bp - last_bp) > radius:
            peaks += 1
            last_chr = row.chr
            last_bp = row.bp

    return peaks
