"""Figure 9: literature-mined external evaluation."""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve

from config import FIGURES, RESULTS

plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight", "font.family": "DejaVu Sans",
                     "axes.titlesize": 7.5})


def main() -> None:
    base = RESULTS / "literature"
    pred = base / "literature_predictions.csv"
    if not pred.exists():
        print("literature predictions missing")
        return
    df = pd.read_csv(pred)
    met = json.loads((base / "metrics.json").read_text()) if (base / "metrics.json").exists() else {}
    nov = df[(df.label == 1) & ~df.in_uniprot]
    neg = df[df.label == 0]

    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6), constrained_layout=True)
    if len(nov) and len(neg):
        y = np.r_[np.ones(len(nov)), np.zeros(len(neg))]
        p = np.r_[nov.p.values, neg.p.values]
        fpr, tpr, _ = roc_curve(y, p)
        ax[0].plot(fpr, tpr, color="#c44e52", lw=1.2,
                   label=f"novel sites (n={len(nov)}), AUC {met.get('novel_vs_neg_auc', float('nan')):.3f}")
        pw = np.r_[nov.pwm.values, neg.pwm.values]
        fpr2, tpr2, _ = roc_curve(y, pw)
        ax[0].plot(fpr2, tpr2, color="#8c8c8c", lw=1.0, ls="--",
                   label=f"PWM sanity check, AUC {met.get('pwm_auc_novel', float('nan')):.3f}")
        ax[0].plot([0, 1], [0, 1], ":", color="k", lw=0.6)
        ax[0].legend(frameon=False, fontsize=5.5, loc="lower right")
    ax[0].set_xlabel("false positive rate")
    ax[0].set_ylabel("true positive rate")
    ax[0].set_title("a  literature sites vs same-protein cysteines", loc="left", fontweight="bold")

    labels = ["all literature\nsites", "already in\nUniProt 2026_03", "novel\n(not annotated)"]
    vals = [met.get("n_literature_sites", len(df[df.label == 1])),
            int((df.label == 1).sum() - len(nov)), len(nov)]
    ax[1].bar(labels, vals, color=["#4c72b0", "#8c8c8c", "#c44e52"], width=0.6)
    for i, v in enumerate(vals):
        ax[1].text(i, v + max(vals) * 0.02, str(v), ha="center", fontsize=6)
    ax[1].set_ylabel("sites")
    ax[1].set_title("b  independence from the training table", loc="left", fontweight="bold")
    fig.savefig(FIGURES / "fig9_literature.png")
    plt.close(fig)
    print("saved fig9")


if __name__ == "__main__":
    main()
