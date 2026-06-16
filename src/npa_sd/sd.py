"""Segmental-duplication track: load the UCSC GRCh38 GenomicSuperDup table and
answer overlap queries quickly.

The index maps each chromosome to two numpy arrays (interval starts and ends),
which lets ``count_overlaps`` test thousands of windows without a Python loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SDIndex = dict[str, tuple[np.ndarray, np.ndarray]]


def load_index(path: str) -> SDIndex:
    """Build an overlap index from the GenomicSuperDup track file.

    Locates the chromosome / start / end columns by header name, falling back to
    the UCSC column order (chrom, chromStart, chromEnd) if the header is absent.
    """
    table = pd.read_csv(path, sep="\t", dtype=str)
    cols = {str(c).lower().lstrip("#"): c for c in table.columns}

    chrom_col = cols.get("chrom") or cols.get("chr")
    start_col = cols.get("chromstart") or cols.get("start")
    end_col = cols.get("chromend") or cols.get("end")

    if not (chrom_col and start_col and end_col):
        # Positional fallback: first column whose values look like "chrN".
        chrom_col = next(
            c for c in table.columns
            if table[c].astype(str).str.startswith("chr").mean() > 0.8
        )
        pos = list(table.columns).index(chrom_col)
        start_col, end_col = table.columns[pos + 1], table.columns[pos + 2]

    table = table[[chrom_col, start_col, end_col]].copy()
    table[start_col] = table[start_col].astype(int)
    table[end_col] = table[end_col].astype(int)

    index: SDIndex = {}
    for chrom, sub in table.groupby(chrom_col):
        index[str(chrom)] = (sub[start_col].to_numpy(), sub[end_col].to_numpy())
    return index


def count_overlaps(index: SDIndex, chrom: str, start: int, end: int) -> int:
    """Number of SD intervals on ``chrom`` overlapping [start, end)."""
    arr = index.get(chrom)
    if arr is None:
        return 0
    starts, ends = arr
    return int(np.count_nonzero((starts < end) & (ends > start)))


def is_sd_window(index: SDIndex, chrom: str, start: int, end: int,
                 tolerance: int = 2) -> bool:
    """Whether a window overlaps more SD intervals than ``tolerance`` allows.

    Windows that fail this test are the ones the SD filter removes.
    """
    return count_overlaps(index, chrom, start, end) > tolerance
