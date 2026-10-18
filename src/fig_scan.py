"""Figure 9: proteome-wide scan and FDR-controlled shortlist."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import FIGURES, RESULTS

plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight", "font.family": "DejaVu Sans",
                     "axes.titlesize": 7.5})


def main() -> None:
    scan = RESULTS / "scan"
    tgt = pd.read_csv(scan / "all_scores.csv")
    dec = pd.read_csv(scan / "decoys_shuffled.csv")
    tab = pd.read_csv(scan / "fdr_summary.csv")

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), constrained_layout=True)
    bins = np.linspace(0, 1, 41)
    ax[0].hist(tgt.p, bins=bins, color="#4c72b0", alpha=0.85, label="human cysteines")
    ax[0].hist(dec.p, bins=bins, color="#c44e52", alpha=0.7, label="shuffled-flank decoys")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("PalmDL score")
    ax[0].set_ylabel("sites (log scale)")
    ax[0].legend(frameon=False, fontsize=6)
    ax[0].set_title("a  score distributions", loc="left", fontweight="bold")

    ax[1].plot(tab.n_predictions, tab.fdr * 100, "o-", color="#c44e52", ms=4)
    for _, r in tab.iterrows():
        ax[1].annotate(f"{r.fdr:.0%}", (r.n_predictions, r.fdr * 100),
                       textcoords="offset points", xytext=(6, -10), fontsize=5.5)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("shortlist size")
    ax[1].set_ylabel("estimated FDR (%)")
    ax[1].set_title("b  target–decoy FDR", loc="left", fontweight="bold")

    x = np.arange(len(tab))
    ax[2].bar(x, tab.known_recovered, color="#8c8c8c", label="known (UniProt)")
    ax[2].bar(x, tab.novel, bottom=tab.known_recovered, color="#55a868", label="novel")
    ax[2].set_xticks(x)
    ax[2].set_xticklabels([f"{int(f * 100)}%" for f in tab.fdr])
    ax[2].set_xlabel("FDR tier")
    ax[2].set_ylabel("sites in shortlist")
    ax[2].legend(frameon=False, fontsize=6)
    ax[2].set_title("c  composition of the shortlist", loc="left", fontweight="bold")
    fig.savefig(FIGURES / "fig9_proteome_scan.png")
    plt.close(fig)
    print("saved fig9")


if __name__ == "__main__":
    main()
