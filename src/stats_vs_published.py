"""Bootstrap comparison of PalmDL's human AUC against published values.

The human cross-validation folds each carry their own model and calibration, so
the comparison uses the mean of the five per-fold human AUCs and a stratified
bootstrap (resampling within folds), matching stats_headline.py.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import RESULTS

OUT = RESULTS / "stats_vs_published.json"


def main() -> None:
    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    d = d[(d.model == "ensemble") & (d.set == "test") & d.is_human.astype(bool)]
    folds = {k: (g.y.values, g.p.values) for k, g in d.groupby("fold")}
    auc = float(np.mean([roc_auc_score(y, p) for y, p in folds.values()]))
    pooled = float(roc_auc_score(d.y.values, d.p.values))

    rng = np.random.default_rng(0)
    vals = []
    for _ in range(4000):
        accs = []
        for y, p in folds.values():
            i = rng.integers(0, len(y), len(y))
            if y[i].min() == y[i].max():
                continue
            accs.append(roc_auc_score(y[i], p[i]))
        if accs:
            vals.append(float(np.mean(accs)))
    vals = np.array(vals)
    out = {
        "human_auc_fold_mean": auc, "human_auc_pooled": pooled,
        "n": int(len(d)), "n_pos": int(d.y.sum()), "n_folds": len(folds),
        "ci95": [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))],
        "p_leq_0900": float((vals <= 0.900).mean()),
        "p_leq_0855": float((vals <= 0.855).mean()),
        "n_boot": int(len(vals)),
        "interpretation": "per-fold mean human AUC with within-fold bootstrap; "
                          "compare with the published values 0.900 (human) and 0.855 (general)",
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
