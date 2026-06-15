#!/usr/bin/env python3
"""Genome-wide test: compare NPA windows in SD-rich vs SD-free bins.

For each chromosome the script downloads the VCF (resumable, cached), classifies
500-kb bins, samples SD-rich and control bins, and scans each. Per-bin results
are written incrementally so an interrupted run resumes where it stopped.

    python scripts/02_genome_wide.py --chroms chr1 chr2 chr7 chr16 --out results
"""

import argparse
import os

import pandas as pd

from npa_sd import config, data, sd, npa, genome_wide, figures


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chroms", nargs="+", default=["chr1", "chr2", "chr7", "chr16"])
    ap.add_argument("--n-each", type=int, default=25, help="bins per class per chromosome")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cache", default="data")
    ap.add_argument("--out", default="results")
    ap.add_argument("--keep-vcf", action="store_true", help="do not delete VCFs after use")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ped = os.path.join(args.cache, "pedigree.txt")
    data.download(config.PEDIGREE_URL, ped)
    trios = data.load_trios(ped)

    sd_track = os.path.join(args.cache, "GRCh38GenomicSuperDup_sort.tab")
    data.download(config.SD_TRACK_URL, sd_track)
    sd_index = sd.load_index(sd_track)

    all_bins = []
    for chrom in args.chroms:
        out_csv = os.path.join(args.out, f"{chrom}_bins.csv")
        done = pd.read_csv(out_csv).to_dict("records") if os.path.exists(out_csv) else []
        completed = {(r["start"], r["end"]) for r in done}

        sd_rich, control = genome_wide.classify_bins(sd_index, chrom)
        chosen = genome_wide.sample_bins(sd_rich, control,
                                         n_each=args.n_each, seed=args.seed)
        pending = [c for c in chosen if (c[0], c[1]) not in completed]
        if not pending:
            print(f"{chrom}: already complete ({len(done)} bins)")
            all_bins.extend(done)
            continue

        print(f"{chrom}: {len(done)}/{len(chosen)} bins done, {len(pending)} to go")
        vcf = data.fetch_vcf(chrom, args.cache)
        if vcf is None:
            print(f"{chrom}: download incomplete, rerun to resume")
            continue

        rows = list(done)
        for start, end, n_sds, label in pending:
            region = os.path.join(args.cache, "_bin.vcf.gz")
            data.extract_region(vcf, chrom, start, end, region)
            windows = npa.scan(region, trios, sd_index=sd_index, filter_sds=False)
            n_unf = len(windows)
            # In a control bin (no SDs) the filter removes nothing.
            n_filt = int((~windows.apply(
                lambda r: sd.is_sd_window(sd_index, r["chrom"], r["start"], r["end"]),
                axis=1)).sum()) if (label == "SD" and n_unf) else n_unf
            rows.append({"chrom": chrom, "start": start, "end": end, "type": label,
                         "n_sds": n_sds, "n_npa_unfiltered": n_unf,
                         "n_npa_filtered": n_filt})
            pd.DataFrame(rows).to_csv(out_csv, index=False)   # checkpoint per bin
            print(f"  {label:7s} {chrom}:{start}-{end}  SDs={n_sds:3d}  "
                  f"NPA {n_unf} -> {n_filt}")

        all_bins.extend(rows)
        if not args.keep_vcf:
            for path in (vcf, vcf + ".tbi"):
                if os.path.exists(path):
                    os.remove(path)

    if all_bins:
        combined = pd.DataFrame(all_bins)
        combined.to_csv(os.path.join(args.out, "sd_control_all.csv"), index=False)
        figures.sd_control_figure(combined, os.path.join(args.out, "fig2_genome_wide.png"))
        summary = combined.groupby("type")[["n_sds", "n_npa_unfiltered", "n_npa_filtered"]].mean()
        print("\nMean per bin by class:")
        print(summary.round(2).to_string())


if __name__ == "__main__":
    main()
