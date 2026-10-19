"""Headline statistics with stratified (per-fold) bootstrap confidence intervals.

Each cross-validation fold is a separate model with its own score calibration, so
pooling raw out-of-fold probabilities across folds mixes scales and biases AUC
downward. We therefore compute every metric inside each fold and average across
folds, and bootstrap within folds to obtain the confidence interval.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, matthews_corrcoef,
                             roc_auc_score)

from config import DATA_PROC, RESULTS
from metrics import find_threshold


def fold_data(te: pd.DataFrame, va: pd.DataFrame, human: bool) -> dict:
    out = {}
    for k in sorted(te.fold.unique()):
        t = te[(te.fold == k)]
        v = va[(va.fold == k)]
        if human:
            t = t[t.is_human.astype(bool)]
            v = v[v.is_human.astype(bool)]
        if len(t) == 0:
            continue
        thr = find_threshold(v.y.values, v.p.values) if len(v) and v.y.nunique() > 1 else 0.5
        out[k] = (t.y.values, t.p.values, (t.p.values >= thr).astype(int))
    return out


def metric(y, p, pred, name: str) -> float:
    if name == "auc":
        return roc_auc_score(y, p)
    if name == "aupr":
        return average_precision_score(y, p)
    return matthews_corrcoef(y, pred)


def point_and_ci(data: dict, name: str, n: int = 2000, seed: int = 0):
    rng = np.random.default_rng(seed)
    if name == "mcc":
        # MCC is a count-based statistic: pool the folds (thresholds stay per fold)
        y_all = np.concatenate([v[0] for v in data.values()])
        pred_all = np.concatenate([v[2] for v in data.values()])
        point = matthews_corrcoef(y_all, pred_all)
        vals = []
        for _ in range(n):
            ys, ps = [], []
            for y, _, pred in data.values():
                i = rng.integers(0, len(y), len(y))
                ys.append(y[i]); ps.append(pred[i])
            vals.append(matthews_corrcoef(np.concatenate(ys), np.concatenate(ps)))
        return float(point), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
    per_fold = [metric(y, p, pred, name) for y, p, pred in data.values()]
    vals = []
    for _ in range(n):
        accs = []
        for y, p, pred in data.values():
            i = rng.integers(0, len(y), len(y))
            if y[i].min() == y[i].max():
                continue
            accs.append(metric(y[i], p[i], pred[i], name))
        if accs:
            vals.append(float(np.mean(accs)))
    return float(np.mean(per_fold)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def delta_auc(data_a: dict, data_b: dict, n: int = 2000, seed: int = 1):
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n):
        accs = []
        for k in data_a:
            if k not in data_b:
                continue
            ya, pa, _ = data_a[k]
            yb, pb, _ = data_b[k]
            assert len(ya) == len(yb)
            i = rng.integers(0, len(ya), len(ya))
            if ya[i].min() == ya[i].max():
                continue
            accs.append(roc_auc_score(ya[i], pa[i]) - roc_auc_score(yb[i], pb[i]))
        if accs:
            diffs.append(float(np.mean(accs)))
    return float(np.mean(diffs)), [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))]


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
        "n_proteins_with_site": int(pos.acc.nunique()),
        "n_organisms": int(pos.organism.nunique()),
    }
    print(json.dumps(stats, indent=2))

    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    out = {"dataset": stats}
    for label, human in (("cv_all", False), ("cv_human", True)):
        te, va = d[d.set == "test"], d[d.set == "val"]
        for m in ["esm", "rf", "cnn", "ensemble"]:
            data = fold_data(te[te.model == m], va[va.model == m], human)
            if not data:
                continue
            out[f"{label}_{m}"] = {name: point_and_ci(data, name)
                                   for name in ("auc", "aupr", "mcc")}
            out[f"{label}_{m}"]["n"] = int(sum(len(v[0]) for v in data.values()))
            out[f"{label}_{m}"]["n_pos"] = int(sum(int(v[0].sum()) for v in data.values()))
            print(label, m, {k: tuple(round(x, 3) for x in v)
                             for k, v in out[f"{label}_{m}"].items() if k in ("auc", "aupr", "mcc")})

    te, va = d[d.set == "test"], d[d.set == "val"]
    ens = fold_data(te[te.model == "ensemble"], va[va.model == "ensemble"], False)
    for other in ("esm", "rf", "cnn"):
        oth = fold_data(te[te.model == other], va[va.model == other], False)
        mean, ci = delta_auc(ens, oth)
        out[f"delta_auc_ensemble_vs_{other}"] = {"mean": mean, "ci": ci}
        print(other, "delta AUC", round(mean, 4), ci)

    (RESULTS / "headline_stats.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("saved", RESULTS / "headline_stats.json")


if __name__ == "__main__":
    main()
