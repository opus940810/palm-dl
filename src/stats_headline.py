"""Headline statistics with bootstrap confidence intervals + dataset summary."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from config import DATA_PROC, RESULTS
from metrics import find_threshold


def per_fold_thresholds(va: pd.DataFrame) -> dict:
    return {k: find_threshold(p.y.values, p.p.values) for k, p in va.groupby("fold")}


def preds_with_thresholds(te: pd.DataFrame, thr: dict) -> np.ndarray:
    out = np.zeros(len(te), dtype=int)
    for k, part in te.groupby("fold"):
        out[part.index.values] = (part.p.values >= thr.get(k, 0.5)).astype(int)
    return out


def boot(y, p, pred, metric, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max():
            continue
        if metric == "auc":
            vals.append(roc_auc_score(y[i], p[i]))
        elif metric == "mcc":
            vals.append(matthews_corrcoef(y[i], pred[i]))
        elif metric == "aupr":
            from sklearn.metrics import average_precision_score
            vals.append(average_precision_score(y[i], p[i]))
    return float(np.mean(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    pos = sites[sites.label == 1]
    stats = {
        "n_sites": int(len(sites)), "n_positive": int(len(pos)),
        "n_positive_experimental": int((pos.source == "uniprot_exp").sum()),
        "n_positive_inferred": int((pos.source == "uniprot_inferred").sum()),
        "n_human_positive": int(pos.is_human.sum()),
        "n_human_positive_proteins": int(pos[pos.is_human].acc.nunique()),
        "n_hard_negative": int((sites.source == "same_protein_cys").sum()),
        "n_background_negative": int((sites.source == "background_cys").sum()),
        "n_unique_proteins": int(sites.acc.nunique()),
        "n_proteins_with_site": int(sites[sites.label == 1].acc.nunique()),
        "n_organisms": int(pos.organism.nunique()),
    }
    print(json.dumps(stats, indent=2))

    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    out = {}
    for label, subset in (("cv_all", "all"), ("cv_human", "human")):
        te, va = d[d.set == "test"], d[d.set == "val"]
        if subset == "human":
            te, va = te[te.is_human.astype(bool)], va[va.is_human.astype(bool)]
        thr = per_fold_thresholds(va)
        for m in ["esm", "rf", "ensemble"]:
            sub = te[te.model == m].reset_index(drop=True)
            vsub = va[va.model == m].reset_index(drop=True)
            t = per_fold_thresholds(vsub)
            pred = preds_with_thresholds(sub, t)
            y, p = sub.y.values, sub.p.values
            out[f"{label}_{m}"] = {
                "auc": boot(y, p, pred, "auc"), "aupr": boot(y, p, pred, "aupr"),
                "mcc": boot(y, p, pred, "mcc"), "n": len(y), "n_pos": int(y.sum()),
            }
            print(label, m, {k: tuple(round(x, 3) for x in v) for k, v in out[f"{label}_{m}"].items()
                             if k in ("auc", "aupr", "mcc")})

    # paired bootstrap: ensemble vs single best model (AUC difference)
    te, va = d[d.set == "test"], d[d.set == "val"]
    ens = te[te.model == "ensemble"].reset_index(drop=True)
    for other in ("esm", "rf"):
        oth = te[te.model == other].reset_index(drop=True)
        assert (ens.y.values == oth.y.values).all() and len(ens) == len(oth)
        rng = np.random.default_rng(1)
        diffs = []
        for _ in range(2000):
            i = rng.integers(0, len(ens), len(ens))
            if ens.y.values[i].min() == ens.y.values[i].max():
                continue
            diffs.append(roc_auc_score(ens.y.values[i], ens.p.values[i])
                         - roc_auc_score(oth.y.values[i], oth.p.values[i]))
        out[f"delta_auc_ensemble_vs_{other}"] = {
            "mean": float(np.mean(diffs)), "ci": [float(np.percentile(diffs, 2.5)),
                                                  float(np.percentile(diffs, 97.5))]}
        print(other, "delta AUC", round(np.mean(diffs), 4), out[f"delta_auc_ensemble_vs_{other}"]["ci"])

    (RESULTS / "headline_stats.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("saved", RESULTS / "headline_stats.json")


if __name__ == "__main__":
    main()
