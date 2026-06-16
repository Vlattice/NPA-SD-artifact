"""Non-parental-allele (NPA) detection in trios.

A trio carries a non-parental allele at a biallelic site when the child genotype
contains an allele present in neither parent. The genome is scanned in fixed-size
sliding windows of SNVs; a window is reported for a trio when its NPA count
reaches ``threshold``.

The per-variant NPA test is vectorized across all trios with numpy. For a
biallelic SNV the test reduces to: the child carries allele a (0 or 1) that is
absent from the pooled parental alleles. This matches the set-based definition
exactly (see tests/test_npa.py) while running fast enough to scan whole
chromosomes.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from cyvcf2 import VCF

from . import config
from .sd import SDIndex, count_overlaps, is_sd_window

COLUMNS = ["child", "pop", "chrom", "start", "end", "npa_count", "n_sds"]


def _trio_indices(samples, trios):
    index = {s: i for i, s in enumerate(samples)}
    valid = [t for t in trios
             if t["child"] in index and t["father"] in index and t["mother"] in index]
    ci = np.array([index[t["child"]] for t in valid])
    fi = np.array([index[t["father"]] for t in valid])
    mi = np.array([index[t["mother"]] for t in valid])
    return valid, ci, fi, mi


def scan(vcf_path: str, trios: list[dict], *,
         window_size: int = config.WINDOW_SIZE,
         step: int = config.STEP,
         threshold: int = config.NPA_THRESHOLD,
         sd_index: Optional[SDIndex] = None,
         filter_sds: bool = False,
         progress_every: int = 0) -> pd.DataFrame:
    """Scan ``vcf_path`` for windows enriched in non-parental alleles.

    Parameters
    ----------
    trios : list of {"child","father","mother","pop"} dicts.
    filter_sds : if True (requires ``sd_index``), drop windows that overlap more
        than ``config.SD_WINDOW_TOLERANCE`` SD intervals.

    Returns one row per (trio, positive window).
    """
    vcf = VCF(vcf_path)

    # Restrict parsing to trio members; cheaper and lower memory.
    members = {s for t in trios for s in (t["child"], t["father"], t["mother"])}
    keep = [s for s in vcf.samples if s in members]
    if keep and len(keep) < len(vcf.samples):
        vcf.set_samples(keep)

    valid, ci, fi, mi = _trio_indices(vcf.samples, trios)
    if not valid:
        return pd.DataFrame(columns=COLUMNS)
    n_trios = len(valid)

    ring = np.zeros((window_size, n_trios), dtype=bool)   # circular NPA buffer
    pos_ring = np.zeros(window_size, dtype=np.int64)
    rows: list[dict] = []
    seen = 0          # variants kept
    chrom = None

    for variant in vcf:
        if not variant.is_snp or len(variant.ALT) != 1:
            continue

        gt = np.asarray(variant.genotypes, dtype=np.int16)[:, :2]
        child, father, mother = gt[ci], gt[fi], gt[mi]
        par0 = (father == 0).any(1) | (mother == 0).any(1)
        par1 = (father == 1).any(1) | (mother == 1).any(1)
        npa = (((child == 0).any(1) & ~par0)
               | ((child == 1).any(1) & ~par1))

        slot = seen % window_size
        ring[slot] = npa
        pos_ring[slot] = variant.POS
        chrom = variant.CHROM
        seen += 1

        if seen >= window_size and (seen - window_size) % step == 0:
            w_start = int(pos_ring[seen % window_size])
            w_end = int(pos_ring[(seen - 1) % window_size])
            if (filter_sds and sd_index is not None
                    and is_sd_window(sd_index, chrom, w_start, w_end,
                                     config.SD_WINDOW_TOLERANCE)):
                continue
            counts = ring.sum(axis=0)
            hits = np.flatnonzero(counts >= threshold)
            if hits.size:
                n_sds = (count_overlaps(sd_index, chrom, w_start, w_end)
                         if sd_index is not None else 0)
                for j in hits:
                    trio = valid[j]
                    rows.append({
                        "child": trio["child"], "pop": trio.get("pop", ""),
                        "chrom": chrom, "start": w_start, "end": w_end,
                        "npa_count": int(counts[j]), "n_sds": n_sds,
                    })

        if progress_every and seen % progress_every == 0:
            print(f"  {seen:,} variants")

    return pd.DataFrame(rows, columns=COLUMNS)
