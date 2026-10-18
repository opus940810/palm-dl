"""Choose ensemble members/weights on validation data only.

Candidates are ranked by pooled validation AUC (never by test performance);
test metrics are reported afterwards for the selected configuration.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from config import RESULTS, SEED
from metrics import find_threshold, metrics_at
from run_cv import balanced_subset

MEMBERS = ["pwm", "lr", "rf", "cnn", "esm"]


def rank01(x: np.ndarray) -> np.ndarray:
    return (rankdata(x) - 0.5) / len(x)


def load(tag_prefix: str) -> dict[str, pd.DataFrame]:
    out = {}
    for m in MEMBERS:
        tag = f"{tag_prefix}_base" if m in ("pwm", "lr", "rf") else f"{tag_prefix}_{m}"
        path = RESULTS / f"preds_{tag}.csv"
        if path.exists():
            d = pd.read_csv(path)
            out[m] = d[d.model == m]
    return out


def eval_weights(models: dict[str, pd.DataFrame], w: dict[str, float]) -> dict:
    aucs, mccs, bal_mccs = [], [], []
    pooled_y, pooled_p = [], []
    val_y, val_p = [], []
    for k in sorted(next(iter(models.values())).fold.unique()):
        p_val, p_test, y_val, y_test = {}, {}, None, None
        for m in models:
            sub = models[m][models[m].fold == k]
            p_val[m] = sub[sub.set == "val"].p.values
            p_test[m] = sub[sub.set == "test"].p.values
            y_val = sub[sub.set == "val"].y.values
            y_test = sub[sub.set == "test"].y.values
        s_val = sum(w[m] * rank01(p_val[m]) for m in models)
        s_test = sum(w[m] * rank01(p_test[m]) for m in models)
        val_y.append(y_val)
        val_p.append(s_val)
        t = find_threshold(y_val, s_val)
        bval = balanced_subset(pd.DataFrame({"label": y_val}), np.zeros(len(y_val)), SEED + k)
        btest = balanced_subset(pd.DataFrame({"label": y_test}), np.zeros(len(y_test)), SEED + 100 + k)
        tb = find_threshold(y_val[bval], s_val[bval])
        aucs.append(metrics_at(y_test, s_test, t)["auc"])
        mccs.append(metrics_at(y_test, s_test, t)["mcc"])
        bal_mccs.append(metrics_at(y_test[btest], s_test[btest], tb)["mcc"])
        pooled_y.append(y_test)
        pooled_p.append(s_test)
    y, p = np.concatenate(pooled_y), np.concatenate(pooled_p)
    vy, vp = np.concatenate(val_y), np.concatenate(val_p)
    return {"val_auc": float(metrics_at(vy, vp, 0.5)["auc"]), "test_auc": float(np.mean(aucs)),
            "mcc": float(np.mean(mccs)), "mcc_sd": float(np.std(mccs)),
            "bal_mcc": float(np.mean(bal_mccs)), "pooled_n": len(y)}


def main() -> None:
    for prefix, label in (("cvh", "human-only CV"), ("cv", "all-species CV")):
        models = load(prefix)
        if len(models) < 2:
            print(prefix, "insufficient members")
            continue
        print(f"\n===== {label} ({prefix}) =====")
        cands = {}
        for r in range(1, len(models) + 1):
            for combo in itertools.combinations(models, r):
                cands["+".join(combo)] = {m: 1.0 / len(combo) if m in combo else 0.0 for m in models}
        rows = []
        for name, w in cands.items():
            r = eval_weights(models, w)
            rows.append({"combo": name, **r})
        t = pd.DataFrame(rows).sort_values("val_auc", ascending=False)
        t.to_csv(RESULTS / f"ensemble_search_{prefix}.csv", index=False)
        print(t.head(12).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
