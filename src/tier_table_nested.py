"""Recompute Table 6 strata and calibration from the nested ensemble predictions.

Reports both pooled out-of-fold values (used in Table 6) and unweighted fold
means, so the manuscript can state which estimator each number uses.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from config import DATA_PROC, RESULTS


def ece(y, p, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, bins - 1)
    total = 0.0
    for b in range(bins):
        m = idx == b
        if m.sum() == 0:
            continue
        total += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(total)


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    meta = sites[["acc", "pos", "evidence", "compartment", "source"]]
    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    d = d[(d.model == "ensemble") & (d.set == "test")].merge(meta, on=["acc", "pos"], how="left")
    exp = d.evidence.astype(str).str.contains("ECO:0000269|ECO:0007744|ECO:0007005|ECO:0006056")

    rows = []
    neg = d[d.y == 0]
    strata = {
        # evidence tiers: positives of the tier plus all negatives
        "experimental": pd.concat([d[exp & (d.y == 1)], neg]),
        "inferred": pd.concat([d[(~exp) & (d.y == 1)], neg]),
    }
    for comp in ["cytoplasmic", "extracellular", "lumenal", "unknown"]:
        strata[comp] = d[d.compartment == comp]

    for name, sub in strata.items():
        if sub.y.nunique() < 2:
            continue
        fold_auc = [roc_auc_score(g.y, g.p) for _, g in sub.groupby("fold") if g.y.nunique() > 1]
        fold_ap = [average_precision_score(g.y, g.p) for _, g in sub.groupby("fold") if g.y.nunique() > 1]
        rows.append({
            "stratum": name, "n_pos": int(sub.y.sum()), "n_rows": len(sub),
            "auc_pooled": float(roc_auc_score(sub.y, sub.p)),
            "ap_pooled": float(average_precision_score(sub.y, sub.p)),
            "auc_foldmean": float(np.mean(fold_auc)), "ap_foldmean": float(np.mean(fold_ap)),
            "n_folds_auc": len(fold_auc),
        })
    tab = pd.DataFrame(rows)
    tab.to_csv(RESULTS / "error_analysis" / "table6_nested.csv", index=False)
    print(tab.round(4).to_string(index=False))

    fold_ece = [ece(g.y.values, g.p.values) for _, g in d.groupby("fold")]
    out = {"ece_foldmean": float(np.mean(fold_ece)), "ece_pooled": ece(d.y.values, d.p.values),
           "n_folds": len(fold_ece)}
    (RESULTS / "error_analysis" / "ece_nested.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
