"""Ensemble out-of-fold predictions and compare models on identical rows.

Each member's scores are converted to within-fold ranks (so that models with
very different score scales can be averaged), the equal-weight average is
thresholded on the inner validation set, and the resulting ensemble is
evaluated on the outer test fold. Row order is shared by all members because
every member is produced by the same fold-splitting code.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from config import DATA_PROC, RESULTS, SEED
from metrics import find_threshold, metrics_at
from run_cv import balanced_subset


def rank01(x: np.ndarray) -> np.ndarray:
    return (rankdata(x) - 0.5) / len(x)


def load(members: list[str], tags: dict[str, str]) -> pd.DataFrame:
    frames = []
    for m in members:
        d = pd.read_csv(RESULTS / f"preds_{tags[m]}.csv")
        frames.append(d[d.model == m])
    return pd.concat(frames)


def attach_human(df: pd.DataFrame) -> pd.DataFrame:
    if "is_human" in df.columns and df.is_human.notna().all():
        return df
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    mp = sites.set_index(["acc", "pos"]).is_human
    df = df.copy()
    df["is_human"] = [bool(mp.get((a, p), False)) for a, p in zip(df.acc, df.pos)]
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", default="pwm,lr,rf,cnn,esm")
    ap.add_argument("--tags", default="", help="model:tag overrides")
    ap.add_argument("--weights", default="", help="comma list of floats, default equal")
    ap.add_argument("--out", default="cv_ensemble")
    args = ap.parse_args()

    members = [m.strip() for m in args.members.split(",")]
    default_tag = "cvh" if args.out.startswith("cvh") else "cv"
    tags = {m: f"{default_tag}_base" if m in ("pwm", "lr", "rf") else f"{default_tag}_{m}" for m in members}
    for item in filter(None, [s.strip() for s in args.tags.split(",")]):
        m, t = item.split(":")
        tags[m] = t
    weights = np.ones(len(members)) / len(members)
    if args.weights:
        weights = np.array([float(w) for w in args.weights.split(",")])
        weights = weights / weights.sum()

    df = attach_human(load(members, tags))
    print(df.groupby(["model", "set"]).size().to_string())
    rows, preds = [], []
    for k in sorted(df.fold.unique()):
        sub = df[df.fold == k]
        p_val = {m: sub[(sub.model == m) & (sub.set == "val")].p.values for m in members}
        p_test = {m: sub[(sub.model == m) & (sub.set == "test")].p.values for m in members}
        ref = sub[(sub.model == members[0]) & (sub.set == "test")]
        y_val = sub[(sub.model == members[0]) & (sub.set == "val")].y.values
        y_test = ref.y.values
        n_val, n_test = len(y_val), len(y_test)
        for m in members:
            assert len(p_val[m]) == n_val and len(p_test[m]) == n_test, f"row mismatch {m} fold {k}"

        def combine(d: dict[str, np.ndarray]) -> np.ndarray:
            return sum(w * rank01(d[m]) for w, m in zip(weights, members))

        s_val, s_test = combine(p_val), combine(p_test)
        # probability calibration on the inner validation fold (Platt scaling)
        from sklearn.linear_model import LogisticRegression
        cal = LogisticRegression(C=1.0, max_iter=1000).fit(s_val.reshape(-1, 1), y_val)
        s_val = cal.predict_proba(s_val.reshape(-1, 1))[:, 1]
        s_test = cal.predict_proba(s_test.reshape(-1, 1))[:, 1]
        t_nat = find_threshold(y_val, s_val)
        bal_val = balanced_subset(pd.DataFrame({"label": y_val}), np.zeros(n_val), SEED + k)
        bal_test = balanced_subset(pd.DataFrame({"label": y_test}), np.zeros(n_test), SEED + 100 + k)
        t_bal = find_threshold(y_val[bal_val], s_val[bal_val])
        m_nat = metrics_at(y_test, s_test, t_nat)
        m_bal = metrics_at(y_test[bal_test], s_test[bal_test], t_bal)
        hard = (y_test == 1) | (ref.source.values == "same_protein_cys")
        m_hard = metrics_at(y_test[hard], s_test[hard], t_nat)
        for tag, m in (("test_natural", m_nat), ("test_balanced", m_bal), ("test_hard", m_hard)):
            rows.append({"model": "ensemble", "fold": k, "protocol": tag,
                         "t_nat": t_nat, "t_bal": t_bal, **m})
        refv = sub[(sub.model == members[0]) & (sub.set == "val")]
        preds.append(pd.DataFrame({"model": "ensemble", "fold": k, "set": "val",
                                   "acc": refv.acc.values, "pos": refv.pos.values,
                                   "y": y_val, "p": s_val, "source": refv.source.values,
                                   "is_human": refv.is_human.astype(bool).values}))
        preds.append(pd.DataFrame({"model": "ensemble", "fold": k, "set": "test",
                                   "acc": ref.acc.values, "pos": ref.pos.values,
                                   "y": y_test, "p": s_test, "source": ref.source.values,
                                   "is_human": ref.is_human.astype(bool).values}))
        print(f"[fold {k}] nat MCC={m_nat['mcc']:.3f} AUC={m_nat['auc']:.3f} | "
              f"bal MCC={m_bal['mcc']:.3f} | hard MCC={m_hard['mcc']:.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / f"metrics_{args.out}.csv", index=False)
    pd.concat(preds).to_csv(RESULTS / f"preds_{args.out}.csv", index=False)
    print(out[["protocol", "mcc", "auc", "aupr", "acc", "sn", "sp", "ece"]].groupby("protocol").mean().round(4))


if __name__ == "__main__":
    main()
