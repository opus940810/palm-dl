"""Additional manuscript figures: cross-species transfer and model interpretation."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import DATA_PROC, FIGURES, RESULTS

plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight", "font.family": "DejaVu Sans",
                     "axes.titlesize": 7.5})


def fig7_matrix() -> None:
    path = RESULTS / "species_matrix_auc.csv"
    if not path.exists():
        print("matrix missing")
        return
    m = pd.read_csv(path, index_col=0)
    order = ["Human", "Mouse", "Rat", "Cattle", "Pig", "A. thaliana", "E. coli", "B. subtilis"]
    order = [s for s in order if s in m.index and s in m.columns]
    m = m.loc[order, order]
    fig, ax = plt.subplots(figsize=(3.9, 3.2), constrained_layout=True)
    im = ax.imshow(m.values, cmap="RdYlBu_r", vmin=0.4, vmax=1.0)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    for i in range(len(order)):
        for j in range(len(order)):
            v = m.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if v < 0.55 or v > 0.93 else "black")
    ax.set_xlabel("test species")
    ax.set_ylabel("training species")
    ax.set_title("AUC of single-species PalmDL models", loc="left", fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=5)
    fig.savefig(FIGURES / "fig7_species_matrix.png")
    plt.close(fig)


def fig8_interpretation() -> None:
    eff_path = RESULTS / "error_analysis" / "mutagenesis_effect.npy"
    fdr_path = RESULTS / "error_analysis" / "fdr_operating_points.csv"
    if not eff_path.exists():
        print("mutagenesis missing")
        return
    eff = np.load(eff_path)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    order = list("IVLFMAWYGPSTCNQDEKRH")
    rows = [aa.index(a) for a in order if a in aa]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.5), constrained_layout=True,
                           gridspec_kw={"width_ratios": [2.1, 1]})
    show = eff[:, 7:24][rows, :]
    im = ax[0].imshow(show, aspect="auto", cmap="RdBu_r", vmin=-0.12, vmax=0.12)
    ax[0].set_yticks(range(len(rows)))
    ax[0].set_yticklabels([aa[i] for i in rows])
    ax[0].set_xticks(range(17))
    ax[0].set_xticklabels(range(-8, 9))
    ax[0].set_xlabel("mutated position relative to cysteine")
    ax[0].set_title("a  in-silico mutagenesis (change in probability)", loc="left", fontweight="bold")
    ax[0].axvline(7, color="k", lw=0.7, ls=":")
    cb = fig.colorbar(im, ax=ax[0], fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=5)

    if fdr_path.exists():
        f = pd.read_csv(fdr_path)
        ax[1].plot(f.test_recall, f.test_precision, "o-", color="#c44e52", ms=3, lw=1)
        for _, r in f.iterrows():
            ax[1].annotate(f"{r.target_fdr:.0%}", (r.test_recall, r.test_precision),
                           textcoords="offset points", xytext=(-9, 5), fontsize=5.5)
        ax[1].set_xlabel("recall of known sites")
        ax[1].set_ylabel("precision")
        ax[1].set_title("b  FDR-controlled operating points", loc="left", fontweight="bold")
        ax[1].set_xlim(0.3, 0.85)
        ax[1].set_ylim(0.6, 1.02)
    fig.savefig(FIGURES / "fig8_interpretation.png")
    plt.close(fig)


if __name__ == "__main__":
    fig7_matrix()
    fig8_interpretation()
    print("saved fig7/fig8")
