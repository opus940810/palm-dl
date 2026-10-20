"""Paired cluster bootstrap of nested-ensemble AUC against each single view."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, RESULTS

rng = np.random.default_rng(20261021)


def main() -> None:
    sp = pd.read_csv(DATA_PROC / "splits_cv.csv")[["acc", "pos", "cluster"]]
    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    d = d[d.set == "test"].merge(sp, on=["acc", "pos"], how="left")
    out = {}
    for other in ("esm", "rf", "cnn"):
        a = d[d.model == "ensemble"].reset_index(drop=True)
        b = d[d.model == other].reset_index(drop=True)
        assert (a.y.values == b.y.values).all() and (a.acc.values == b.acc.values).all()
        diffs = []
        for _ in range(1000):
            vals = []
            for k in sorted(a.fold.unique()):
                m = (a.fold == k).values
                c = a.cluster.values[m]
                codes, uniq = pd.factorize(c)
                mult = rng.multinomial(len(uniq), np.full(len(uniq), 1.0 / len(uniq)))[codes]
                idx = np.repeat(np.arange(m.sum()), mult)
                ya = a.y.values[m][idx]
                if ya.min() == ya.max():
                    continue
                vals.append(roc_auc_score(ya, a.p.values[m][idx]) -
                            roc_auc_score(ya, b.p.values[m][idx]))
            if vals:
                diffs.append(float(np.mean(vals)))
        per_fold = [roc_auc_score(a.y.values[(a.fold == k).values], a.p.values[(a.fold == k).values])
                    - roc_auc_score(b.y.values[(b.fold == k).values], b.p.values[(b.fold == k).values])
                    for k in sorted(a.fold.unique())]
        out[other] = {"mean": float(np.mean(per_fold)),
                      "ci_cluster": [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))]}
        print(other, round(out[other]["mean"], 4), [round(x, 4) for x in out[other]["ci_cluster"]])
    (RESULTS / "delta_nested.json").write_text(json.dumps(out, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
