"""Evaluate an apparent NPA signal that falls outside annotated SDs.

A collapsed paralog (two near-identical copies mapped to one locus) inflates
heterozygosity far above what the allele frequency predicts. Comparing observed
heterozygosity with the frequency-based expectation at each site is a cheap test
for that failure mode, using the data already at hand rather than a re-alignment.

Note: a clean result here rules out collapsed paralogs only. It does not clear a
locus that sits in subtelomeric or otherwise low-mappability sequence, where
short-read genotypes can be unreliable for other reasons.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from cyvcf2 import VCF

HET_EXCESS = 0.55  # observed het above this flags a likely collapsed paralog


def heterozygosity_profile(vcf_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Per-site allele frequency and heterozygous fraction across all samples.

    ``vcf_path`` must be indexed (.tbi/.csi) so the region can be queried.
    """
    vcf = VCF(vcf_path)
    rows = []
    for variant in vcf(f"{chrom}:{start}-{end}"):
        if not variant.is_snp or len(variant.ALT) != 1:
            continue
        gt = np.asarray(variant.genotypes, dtype=np.int16)[:, :2]
        called = (gt >= 0).all(axis=1)
        alleles = gt[called]
        if alleles.size == 0:
            continue
        het = float((alleles[:, 0] != alleles[:, 1]).mean())
        af = float(alleles.sum() / (2 * len(alleles)))
        rows.append({
            "pos": int(variant.POS),
            "af": af,
            "het": het,
            "het_expected": 2 * af * (1 - af),
            "het_excess": het > HET_EXCESS,
        })
    return pd.DataFrame(rows)
