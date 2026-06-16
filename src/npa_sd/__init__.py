"""npa_sd: reanalysis of non-parental-allele clusters in 1000 Genomes trios.

Public API:
    load_trios, fetch_vcf, extract_region   (data)
    load_index, count_overlaps, is_sd_window (sd)
    scan                                     (npa)
    classify_bins, sample_bins               (genome_wide)
    heterozygosity_profile                   (candidate)
"""

from . import config
from .data import load_trios, fetch_vcf, extract_region, download
from .sd import load_index, count_overlaps, is_sd_window
from .npa import scan
from .genome_wide import classify_bins, sample_bins
from .candidate import heterozygosity_profile

__all__ = [
    "config", "load_trios", "fetch_vcf", "extract_region", "download",
    "load_index", "count_overlaps", "is_sd_window", "scan",
    "classify_bins", "sample_bins", "heterozygosity_profile",
]

__version__ = "0.1.0"
