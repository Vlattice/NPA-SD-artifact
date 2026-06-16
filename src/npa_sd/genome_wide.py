"""Genome-wide test: do NPA windows concentrate in segmental-duplication regions?

A chromosome is tiled into fixed bins. Bins with many SD intervals ("SD-rich")
are compared against bins with none ("control"). If the reported NPA signal is an
SD artifact, it should appear in the SD-rich bins and be absent from the controls.
"""

from __future__ import annotations

import random

from . import config
from .sd import SDIndex, count_overlaps


def classify_bins(index: SDIndex, chrom: str, *,
                  bin_size: int = config.SD_BIN_SIZE,
                  sd_rich_min: int = config.SD_RICH_MIN):
    """Split a chromosome into SD-rich and SD-free bins.

    Returns (sd_rich, control), each a list of (start, end, n_sds). Bins with an
    intermediate SD count (between 1 and ``sd_rich_min``) are left out so the two
    groups contrast cleanly.
    """
    length = config.CHROM_LENGTHS[chrom]
    sd_rich, control = [], []
    for start in range(0, length - bin_size, bin_size):
        end = start + bin_size
        n = count_overlaps(index, chrom, start, end)
        if n >= sd_rich_min:
            sd_rich.append((start, end, n))
        elif n == 0:
            control.append((start, end, n))
    return sd_rich, control


def sample_bins(sd_rich, control, *, n_each: int = 25, seed: int = 0):
    """Randomly pick ``n_each`` bins from each group. Returns (start, end, n_sds, label)."""
    rng = random.Random(seed)
    sd_rich = list(sd_rich)
    control = list(control)
    rng.shuffle(sd_rich)
    rng.shuffle(control)
    chosen = [(*b, "SD") for b in sd_rich[:n_each]]
    chosen += [(*b, "control") for b in control[:n_each]]
    return chosen
