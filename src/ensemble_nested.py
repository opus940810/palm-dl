"""Nested ensemble selection: choose members inside each outer fold.

The previous protocol ranked the 31 member subsets by pooled validation AUC
across all folds and then applied one global combination to every test fold.
That lets each fold's validation rows influence other folds' test predictions.
Here the combination is chosen using only the validation rows of the same outer
fold, calibrated on those rows, and evaluated on that fold's test rows.
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from config import RESULTS, SEED
from metrics import find_threshold, metrics_at
from run_cv import balanced_subset

MEMBERS = ["pwm", "lr", "rf", "cnn", "esm"]


def rank01(x: np.ndarray) -> np.ndarray:
    return (rankdata(x) - 0.5) / len(x)


def load(prefix: str) -> dict[str, pd.DataFrame]:
    out = {}
    for m in MEMBERS:
        tag = f"{prefix}_base" if m in ("pwm", "lr", "rf") else f"{prefix}_{m}"
        path = RESULTS / f"preds_{tag}.csv"
        if path.exists():
            d = pd.read_csv(path)
            out[m] = d[d.model == m]
    return out


def combos(models: dict[str, pd.DataFrame]):
    names = list(models)
    for r in range(1, len(names) + 1):
        for c in itertools.combinations(names, r):
            yield c


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="cv")
    ap.add_argument("--out", default="cv_ensemble")
    args = ap.parse_args()

    models = load(args.prefix)
    if len(models) < 2:
        raise SystemExit(f"need at least two members, found {list(models)}")
    print("members:", list(models))

    rows, preds, choices = [], [], []
    for k in sorted(next(iter(models.values())).fold.unique()):
        y_val = None
        y_test = None
        p_val, p_test = {}, {}
        for m, d in models.items():
            sub = d[d.fold == k]
            v = sub[sub.set == "val"]
            t = sub[sub.set == "test"]
            p_val[m] = v.p.values
            p_test[m] = t.p.values
            y_val = v.y.values
            y_test = t.y.values
            ref_v, ref_t = v, t
        # choose the combination on this fold's validation rows only
        best = None
        for combo in combos(models):
            s = sum(rank01(p_val[m]) for m in combo) / len(combo)
            auc = roc_auc_score(y_val, s)
            key = (round(float(auc), 6), -len(combo), combo)
            if best is None or key > best[0]:
                best = (key, combo, float(auc))
        _, combo, val_auc = best
        s_val = sum(rank01(p_val[m]) for m in combo) / len(combo)
        s_test = sum(rank01(p_test[m]) for m in combo) / len(combo)
        cal = LogisticRegression(C=1.0, max_iter=1000).fit(s_val.reshape(-1, 1), y_val)
        q_val = cal.predict_proba(s_val.reshape(-1, 1))[:, 1]
        q_test = cal.predict_proba(s_test.reshape(-1, 1))[:, 1]
        t_nat = find_threshold(y_val, q_val)
        bal_val = balanced_subset(pd.DataFrame({"label": y_val}), np.zeros(len(y_val)), SEED + k)
        bal_test = balanced_subset(pd.DataFrame({"label": y_test}), np.zeros(len(y_test)), SEED + 100 + k)
        t_bal = find_threshold(y_val[bal_val], q_val[bal_val])
        m_nat = metrics_at(y_test, q_test, t_nat)
        m_bal = metrics_at(y_test[bal_test], q_test[bal_test], t_bal)
        rows.append({"model": "ensemble", "fold": k, "protocol": "test_natural", "t_nat": t_nat,
                     "t_bal": t_bal, **m_nat, "members": "+".join(combo), "val_auc": val_auc})
        rows.append({"model": "ensemble", "fold": k, "protocol": "test_balanced", "t_nat": t_nat,
                     "t_bal": t_bal, **m_bal, "members": "+".join(combo), "val_auc": val_auc})
        choices.append({"fold": k, "members": "+".join(combo), "val_auc": val_auc,
                        "test_auc": m_nat["auc"]})
        preds.append(pd.DataFrame({"model": "ensemble", "fold": k, "set": "val",
                                   "acc": ref_v.acc.values, "pos": ref_v.pos.values, "y": y_val,
                                   "p": q_val, "source": ref_v.source.values,
                                   "is_human": ref_v.is_human.astype(bool).values}))
        preds.append(pd.DataFrame({"model": "ensemble", "fold": k, "set": "test",
                                   "acc": ref_t.acc.values, "pos": ref_t.pos.values, "y": y_test,
                                   "p": q_test, "source": ref_t.source.values,
                                   "is_human": ref_t.is_human.astype(bool).values}))
        print(f"fold {k}: {combo} val AUC {val_auc:.4f} test AUC {m_nat['auc']:.4f} "
              f"MCC {m_nat['mcc']:.3f}")

    met = pd.DataFrame(rows)
    met.to_csv(RESULTS / f"metrics_{args.out}.csv", index=False)
    pd.concat(preds).to_csv(RESULTS / f"preds_{args.out}.csv", index=False)
    pd.DataFrame(choices).to_csv(RESULTS / f"ensemble_nested_{args.prefix}.csv", index=False)
    print("\nnested mean test AUC:", round(float(met[met.protocol == "test_natural"].auc.mean()), 5))
    print(met[met.protocol == "test_natural"][["fold", "members", "auc", "mcc"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
