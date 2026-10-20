"""Headline statistics: site-level and cluster-level bootstrap intervals.

Each cross-validation fold is a separate model with its own calibration, so
discrimination is computed inside folds and averaged (unweighted fold mean);
MCC pools confusion counts at per-fold thresholds. Two bootstrap schemes are
reported: resampling sites (the optimistic interval used previously) and
resampling protein clusters (accounts for within-cluster correlation).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, matthews_corrcoef,
                             roc_auc_score)

from config import DATA_PROC, RESULTS
from metrics import find_threshold


def load_clusters() -> dict:
    sp = pd.read_csv(DATA_PROC / "splits_cv.csv")[["acc", "pos", "cluster"]]
    return {(a, p): c for a, p, c in zip(sp.acc, sp.pos, sp.cluster)}


def fold_data(te: pd.DataFrame, va: pd.DataFrame, human: bool) -> dict:
    out = {}
    for k in sorted(te.fold.unique()):
        t = te[te.fold == k]
        v = va[va.fold == k]
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


def site_ci(data: dict, name: str, n: int = 2000, seed: int = 0):
    rng = np.random.default_rng(seed)
    if name == "mcc":
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


def cluster_ci(data: dict, clusters: dict, name: str, n: int = 1000, seed: int = 20261020):
    """Resample protein clusters within each fold (clusters are the sampling unit).

    Uses a multinomial draw of cluster multiplicities and np.repeat, which is
    equivalent to resampling clusters with replacement but avoids per-cluster
    index construction.
    """
    rng = np.random.default_rng(seed)
    packed = {}
    for k, (y, p, pred) in data.items():
        c = np.asarray(clusters[k])
        codes, uniq = pd.factorize(c)
        packed[k] = (y, p, pred, codes, len(uniq))
    if name == "mcc":
        y_all = np.concatenate([v[0] for v in data.values()])
        pred_all = np.concatenate([v[2] for v in data.values()])
        point = matthews_corrcoef(y_all, pred_all)
        vals = []
        for _ in range(n):
            ys, ps = [], []
            for y, _, pred, codes, nclus in packed.values():
                mult = rng.multinomial(nclus, np.full(nclus, 1.0 / nclus))[codes]
                idx = np.repeat(np.arange(len(y)), mult)
                ys.append(y[idx]); ps.append(pred[idx])
            vals.append(matthews_corrcoef(np.concatenate(ys), np.concatenate(ps)))
        return float(point), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
    per_fold = [metric(y, p, pred, name) for y, p, pred in data.values()]
    vals = []
    for _ in range(n):
        accs = []
        for y, p, pred, codes, nclus in packed.values():
            mult = rng.multinomial(nclus, np.full(nclus, 1.0 / nclus))[codes]
            idx = np.repeat(np.arange(len(y)), mult)
            if y[idx].min() == y[idx].max():
                continue
            accs.append(metric(y[idx], p[idx], pred[idx], name))
        if accs:
            vals.append(float(np.mean(accs)))
    return float(np.mean(per_fold)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


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

    cmap = load_clusters()
    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    d["cluster"] = [cmap.get((a, p), f"{a}_{p}") for a, p in zip(d.acc, d.pos)]

    out = {"dataset": stats}
    for label, human in (("cv_all", False), ("cv_human", True)):
        te, va = d[d.set == "test"], d[d.set == "val"]
        for m in ["esm", "rf", "cnn", "ensemble"]:
            data = fold_data(te[te.model == m], va[va.model == m], human)
            if not data:
                continue
            cl = {}
            for k in data:
                mask = (te.model == m) & (te.fold == k)
                if human:
                    mask = mask & te.is_human.astype(bool)
                cl[k] = te[mask].cluster.values
            entry = {}
            for name in ("auc", "aupr", "mcc"):
                entry[name] = site_ci(data, name)
                entry[f"{name}_cluster"] = cluster_ci(data, cl, name)
            entry["n"] = int(sum(len(v[0]) for v in data.values()))
            entry["n_pos"] = int(sum(int(v[0].sum()) for v in data.values()))
            out[f"{label}_{m}"] = entry
            print(label, m, "AUC site", tuple(round(x, 3) for x in entry["auc"]),
                  "cluster", tuple(round(x, 3) for x in entry["auc_cluster"]))

    (RESULTS / "headline_stats.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("saved", RESULTS / "headline_stats.json")


if __name__ == "__main__":
    main()
