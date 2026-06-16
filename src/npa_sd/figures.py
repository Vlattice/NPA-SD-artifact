"""Figures used in the reanalysis, generated from result tables.

These render the two data figures (chromosome-3 reproduction and the genome-wide
SD-vs-control comparison). They take DataFrames produced by the scripts so a
figure can be regenerated without rerunning the scan.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

_RED, _BLUE, _GREY, _INK = "#c0392b", "#2c6e9b", "#9aa5ad", "#1f2a33"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.8, "font.size": 10,
})


def chr3_figure(windows: pd.DataFrame, out_path: str,
                title: str = "chr3:75.2-75.7 Mb") -> str:
    """NPA windows along the chromosome-3 locus, coloured by SD overlap."""
    fig, ax = plt.subplots(figsize=(7, 3.6))
    if len(windows):
        mid = (windows["start"] + windows["end"]) / 2 / 1e6
        sc = ax.scatter(mid, windows["npa_count"], c=windows["n_sds"],
                        cmap="OrRd", s=34, edgecolor=_INK, linewidth=0.3, zorder=3)
        fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.05,
                     label="SD intervals per window")
    ax.axhline(5, ls="--", lw=0.9, color=_GREY)
    ax.set_xlabel("Position on chr3 (Mb)")
    ax.set_ylabel("NPAs per window")
    ax.set_title(f"{title}: {len(windows)} positive windows", loc="left",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def sd_control_figure(bins: pd.DataFrame, out_path: str) -> str:
    """SD count vs NPA windows per bin, plus the per-class mean."""
    sd = bins[bins["type"] == "SD"]
    ctrl = bins[bins["type"] == "control"]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6),
                           gridspec_kw={"width_ratios": [1.4, 1]})
    ax[0].scatter(sd["n_sds"], sd["n_npa_unfiltered"], c=_RED, s=26, alpha=.65,
                  edgecolor="none", label="SD-rich bins")
    ax[0].scatter(ctrl["n_sds"], ctrl["n_npa_unfiltered"], c=_BLUE, s=26,
                  alpha=.75, edgecolor="none", label="control bins")
    ax[0].set_xlabel("SD intervals per 500-kb bin")
    ax[0].set_ylabel("Positive NPA windows")
    ax[0].set_title("a  SD density vs NPA signal", loc="left", fontweight="bold")
    ax[0].legend(frameon=False, fontsize=8)

    means = bins.groupby("type")["n_npa_unfiltered"].mean().reindex(["SD", "control"])
    ax[1].bar([0, 1], [means["SD"], means["control"]], color=[_RED, _BLUE], width=0.6)
    ax[1].set_xticks([0, 1])
    ax[1].set_xticklabels(["SD-rich", "control"])
    ax[1].set_ylabel("Mean NPA windows per bin")
    ax[1].set_title("b  Per-class mean", loc="left", fontweight="bold")
    for i, v in enumerate([means["SD"], means["control"]]):
        ax[1].text(i, v + 0.2, f"{v:.1f}", ha="center", fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path
