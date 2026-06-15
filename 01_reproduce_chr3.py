#!/usr/bin/env python3
"""Reproduce the strongest reported chromosome-3 signal, then apply the SD filter.

Expected outcome: ~80 positive windows at chr3:75.2-75.7 Mb without the filter,
0 after it, and 0 at the SD-free FRA3B locus either way.

    python scripts/01_reproduce_chr3.py --cache data --out results
"""

import argparse
import os

from npa_sd import config, data, sd, npa, figures


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", default="data", help="directory for downloaded inputs")
    ap.add_argument("--out", default="results", help="directory for outputs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("Loading pedigree and SD track...")
    ped = os.path.join(args.cache, "pedigree.txt")
    data.download(config.PEDIGREE_URL, ped)
    trios = data.load_trios(ped)
    print(f"  complete trios: {len(trios)}")

    sd_track = os.path.join(args.cache, "GRCh38GenomicSuperDup_sort.tab")
    data.download(config.SD_TRACK_URL, sd_track)
    sd_index = sd.load_index(sd_track)

    print("Fetching chr3 (resumable)...")
    vcf = data.fetch_vcf("chr3", args.cache)

    for label, (chrom, start, end) in [("hotspot", config.CHR3_SD_HOTSPOT),
                                       ("fra3b", config.CHR3_FRA3B)]:
        region = os.path.join(args.cache, f"chr3_{label}.vcf.gz")
        data.extract_region(vcf, chrom, start, end, region)

        unfiltered = npa.scan(region, trios, sd_index=sd_index, filter_sds=False)
        filtered = npa.scan(region, trios, sd_index=sd_index, filter_sds=True)

        unfiltered.to_csv(os.path.join(args.out, f"chr3_{label}_unfiltered.csv"), index=False)
        filtered.to_csv(os.path.join(args.out, f"chr3_{label}_filtered.csv"), index=False)

        n0, n1 = len(unfiltered), len(filtered)
        drop = 100 * (n0 - n1) / n0 if n0 else 0
        print(f"chr3 {label}: {n0} windows -> {n1} after SD filter ({drop:.0f}% removed)")

        if label == "hotspot":
            figures.chr3_figure(unfiltered,
                                os.path.join(args.out, "fig1_chr3.png"))


if __name__ == "__main__":
    main()
