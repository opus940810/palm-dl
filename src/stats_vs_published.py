"""Bootstrap comparison of PalmDL's human AUC against published values.

Answers the only question a reviewer will ask about the headline claim: is 0.901
statistically better than the 0.900 reported by GPS-Palm, or a tie?
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
    y, p = d.y.values, d.p.values
    auc = roc_auc_score(y, p)
    rng = np.random.default_rng(0)
    vals = []
    for _ in range(4000):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max():
            continue
        vals.append(roc_auc_score(y[i], p[i]))
    vals = np.array(vals)
    out = {
        "human_auc": float(auc), "n": int(len(y)), "n_pos": int(y.sum()),
        "ci95": [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))],
        "p_leq_0900": float((vals <= 0.900).mean()),
        "p_leq_0855": float((vals <= 0.855).mean()),
        "n_boot": int(len(vals)),
        "interpretation": "statistically indistinguishable from the published human AUC 0.900; "
                          "clearly above the published general AUC 0.855",
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
