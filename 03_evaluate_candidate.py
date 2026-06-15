#!/usr/bin/env python3
"""Evaluate the one apparent NPA signal outside annotated SDs (chr16:88.5-89 Mb).

Reports per-site heterozygosity against the frequency-based expectation. A clean
profile rules out a collapsed paralog but does not validate the locus: it sits in
subtelomeric 16q24.3, where short-read genotypes are unreliable, so the signal
cannot be confirmed without long-read data.

    python scripts/03_evaluate_candidate.py --cache data --out results
"""

import argparse
import os

from npa_sd import config, data, candidate


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    chrom, start, end = config.CANDIDATE_LOCUS
    vcf = data.fetch_vcf(chrom, args.cache)

    profile = candidate.heterozygosity_profile(vcf, chrom, start, end)
    profile.to_csv(os.path.join(args.out, "candidate_het_profile.csv"), index=False)

    n_excess = int(profile["het_excess"].sum()) if len(profile) else 0
    print(f"Sites examined: {len(profile)}")
    print(f"Sites with heterozygosity excess (> {candidate.HET_EXCESS:.0%}): {n_excess}")
    if n_excess == 0:
        print("No collapsed-paralog signature. The locus is still subtelomeric "
              "(16q24.3) and unverifiable from short reads; treat as a probable "
              "mapping artifact pending long-read confirmation.")
    else:
        print("Heterozygosity excess present: consistent with a collapsed paralog "
              "(mapping artifact).")


if __name__ == "__main__":
    main()
