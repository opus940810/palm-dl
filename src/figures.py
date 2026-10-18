"""Manuscript figures (300 dpi PNG)."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve

from config import DATA_PROC, FIGURES, RESULTS

plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight", "font.family": "DejaVu Sans",
                     "axes.titlesize": 7.5, "legend.fontsize": 6})
ORG = {"9606": "Human", "10090": "Mouse", "10116": "Rat", "9913": "Cattle", "83333": "E. coli",
       "3702": "A. thaliana", "9823": "Pig", "224308": "B. subtilis", "9031": "Chicken",
       "9615": "Dog", "559292": "S. cerevisiae", "7955": "Zebrafish"}
COL = {"pwm": "#8c8c8c", "lr": "#4c72b0", "rf": "#55a868", "cnn": "#c44e52",
       "esm": "#8172b2", "ensemble": "#000000"}
LAB = {"pwm": "PWM (CSS-Palm-style)", "lr": "Logistic regression", "rf": "Random forest",
       "cnn": "PalmDL-CNN", "esm": "ESM-2 embedding", "ensemble": "PalmDL ensemble"}


def fig1() -> None:
    s = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    pos = s[s.label == 1]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3), constrained_layout=True)
    vc = pos.organism.value_counts()
    head = vc.head(10)
    names = [ORG.get(str(i), "other") for i in head.index]
    others = int(vc.iloc[10:].sum())
    colors = ["#c44e52" if n == "Human" else ("#4c72b0" if n == "Mouse" else "#c9c9c9") for n in names]
    ax[0].barh(range(len(head))[::-1], head.values, color=colors, height=0.75)
    ax[0].set_yticks(range(len(head))[::-1])
    ax[0].set_yticklabels(names)
    ax[0].set_xlabel(f"S-palmitoyl sites\n(other {len(vc) - 10} species: {others:,})")
    ax[0].set_xlim(0, head.max() * 1.18)
    ax[0].set_title("a  species", loc="left", fontweight="bold")

    t = pos.source.value_counts()
    ax[1].bar(["experimental", "inferred"], [t.get("uniprot_exp", 0), t.get("uniprot_inferred", 0)],
              color=["#c44e52", "#c9c9c9"], width=0.6)
    for i, v in enumerate([t.get("uniprot_exp", 0), t.get("uniprot_inferred", 0)]):
        ax[1].text(i, v * 1.03, str(v), ha="center", fontsize=6)
    ax[1].set_yscale("log")
    ax[1].set_ylabel("sites (log scale)")
    ax[1].set_title("b  evidence tier", loc="left", fontweight="bold")

    comp = pos.compartment.value_counts()
    order = [c for c in ["cytoplasmic", "unknown", "extracellular", "lumenal", "other", "signal_peptide"] if c in comp.index]
    ax[2].bar(range(len(order)), [comp[c] for c in order], color="#55a868", width=0.65)
    ax[2].set_xticks(range(len(order)))
    ax[2].set_xticklabels([c.replace("_", " ") for c in order], rotation=35, ha="right")
    ax[2].set_ylabel("sites")
    ax[2].set_title("c  annotated topology", loc="left", fontweight="bold")
    fig.savefig(FIGURES / "fig1_dataset.png")
    plt.close(fig)


def fig2() -> None:
    d = pd.read_csv(RESULTS / "leakage_ladder.csv")
    labels = ["random\n(site-level)", "protein-\nlevel", "cluster-level\n(this work)", "species holdout\n(human test)"]
    x = np.arange(len(d))
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6), constrained_layout=True)
    ax[0].bar(x - 0.19, d.auc, 0.37, color="#4c72b0", label="AUC")
    ax[0].bar(x + 0.19, d.mcc_bal, 0.37, color="#c44e52", label="MCC (balanced)")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(labels)
    ax[0].set_ylim(0, 1.08)
    ax[0].legend(frameon=False, loc="upper right", ncol=1)
    ax[0].set_ylabel("score")
    ax[0].set_title("a  performance depends on the split", loc="left", fontweight="bold")

    ax[1].bar(x, d.leaked_pos_windows * 100, color="#8172b2", width=0.6)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(labels)
    ax[1].set_ylim(0, 55)
    ax[1].set_ylabel("% test positives with an\nidentical training window")
    ax[1].set_title("b  measured leakage", loc="left", fontweight="bold")
    for i, v in enumerate(d.leaked_pos_windows * 100):
        ax[1].text(i, v + 1.2, f"{v:.0f}%", ha="center", fontsize=6)
    fig.savefig(FIGURES / "fig2_leakage.png")
    plt.close(fig)


def fig3() -> None:
    path = RESULTS / "preds_cv_all.csv"
    if not path.exists():
        return
    d = pd.read_csv(path)
    d = d[d.set == "test"]
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.9), constrained_layout=True)
    for m in ["pwm", "lr", "rf", "cnn", "esm", "ensemble"]:
        sub = d[d.model == m]
        if sub.empty:
            continue
        fpr, tpr, _ = roc_curve(sub.y, sub.p)
        ax[0].plot(fpr, tpr, lw=1.1, color=COL[m], label=LAB[m])
        pr, rc, _ = precision_recall_curve(sub.y, sub.p)
        ax[1].plot(rc, pr, lw=1.1, color=COL[m])
    ax[0].plot([0, 1], [0, 1], "--", color="k", lw=0.6)
    ax[0].set_xlabel("false positive rate")
    ax[0].set_ylabel("true positive rate")
    ax[0].legend(frameon=False, loc="lower right")
    ax[0].set_title("a  ROC, cluster-held-out CV", loc="left", fontweight="bold")
    ax[1].set_xlabel("recall")
    ax[1].set_ylabel("precision")
    ax[1].set_title("b  precision-recall", loc="left", fontweight="bold")
    fig.savefig(FIGURES / "fig3_roc.png")
    plt.close(fig)


def fig4() -> None:
    rows = []
    for proto, f, auc_key, mcc_key in (
            ("multi-species\ntraining", "table2_cv_human.csv", "nat_auc", "bal_mcc"),
            ("human-only\ntraining", "table3_cv_human_only.csv", "nat_auc", "bal_mcc"),
            ("species\nholdout", "table4_species.csv", "auc_natural", "mcc_balanced")):
        p = RESULTS / "tables" / f
        if not p.exists():
            continue
        t = pd.read_csv(p)
        for _, r in t.iterrows():
            rows.append({"protocol": proto, "model": r["model"],
                         "auc": r[auc_key], "mcc": r[mcc_key]})
    t = pd.DataFrame(rows)
    protos = list(dict.fromkeys(t.protocol))
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.7), constrained_layout=True)
    for j, metric in enumerate(["auc", "mcc"]):
        for m in ["pwm", "rf", "cnn", "esm", "ensemble"]:
            sub = t[t.model == m]
            if sub.empty:
                continue
            xs = [protos.index(p) for p in sub.protocol]
            ax[j].scatter(xs, sub[metric].values, s=26, color=COL[m], label=LAB[m], zorder=3)
        ax[j].set_xticks(range(len(protos)))
        ax[j].set_xticklabels(protos)
        ax[j].set_ylim(0, 1.0)
        ax[j].grid(axis="y", lw=0.3, alpha=0.5)
        ax[j].set_ylabel("AUC" if metric == "auc" else "MCC (balanced)")
        ax[j].set_title(("a  discrimination on human sites" if metric == "auc"
                         else "b  threshold-dependent accuracy"), loc="left", fontweight="bold")
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, -0.06))
    fig.savefig(FIGURES / "fig4_protocols.png")
    plt.close(fig)


def fig5() -> None:
    path = RESULTS / "preds_cv_all.csv"
    if not path.exists():
        return
    d = pd.read_csv(path)
    d = d[d.set == "test"]
    models = [m for m in ["rf", "cnn", "esm", "ensemble"] if m in set(d.model)]
    if not models:
        return
    fig, ax = plt.subplots(1, len(models), figsize=(1.85 * len(models) + 0.5, 2.05), squeeze=False,
                           constrained_layout=True)
    for i, m in enumerate(models):
        sub = d[d.model == m]
        bins = np.linspace(0, 1, 11)
        idx = np.digitize(sub.p, bins) - 1
        xs, ys = [], []
        for b in range(10):
            mask = idx == b
            if mask.sum() > 20:
                xs.append(sub.p[mask].mean())
                ys.append(sub.y[mask].mean())
        a = ax[0, i]
        a.plot([0, 1], [0, 1], "--", color="k", lw=0.6)
        a.plot(xs, ys, "o-", color=COL.get(m, "#333"), ms=2.5, lw=1)
        a.set_title(LAB.get(m, m), fontsize=6.5)
        a.set_xlabel("predicted")
        if i == 0:
            a.set_ylabel("observed")
    fig.savefig(FIGURES / "fig5_calibration.png")
    plt.close(fig)


def fig6() -> None:
    """Position-specific amino-acid log-odds around the central cysteine."""
    from config import AA

    s = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    pos = s[(s.label == 1) & s.window.str.len().eq(31)].window
    pos = pos[pos.str[15].eq("C")]  # 1% of UniProt annotations do not map to a cysteine
    neg = s[(s.label == 0) & (s.source == "background_cys") & s.window.str.len().eq(31)].window
    # equal sample sizes so that additive pseudocounts cannot create spurious log-odds
    n = min(len(pos), len(neg))
    pos = pos.sample(n=n, random_state=0)
    neg = neg.sample(n=n, random_state=0)
    idx = {a: i for i, a in enumerate(AA)}
    fp = np.full((31, 20), 0.5)
    fn = np.full((31, 20), 0.5)
    for ws, mat in ((pos, fp), (neg, fn)):
        for w in ws:
            for j, ch in enumerate(w):
                if ch in idx:
                    mat[j, idx[ch]] += 1
    fp /= fp.sum(1, keepdims=True)
    fn /= fn.sum(1, keepdims=True)
    lo = np.log2(fp / fn)
    order = list("IVLFMAWYGPSTCNQDEKRH")
    rows = [AA.index(a) for a in order if a in AA]
    fig, ax = plt.subplots(figsize=(4.8, 2.3), constrained_layout=True)
    im = ax.imshow(lo[8:24][:, rows].T, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([AA[i] for i in rows])
    ax.set_xticks(range(16))
    ax.set_xticklabels(range(-7, 9))
    ax.set_xlabel("position relative to cysteine")
    ax.axvline(7, color="k", lw=0.8, ls=":")
    ax.set_title("log2 odds, positive vs background", loc="left", fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.ax.tick_params(labelsize=5)
    fig.savefig(FIGURES / "fig6_motif.png")
    plt.close(fig)


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    print("figures written to", FIGURES)
