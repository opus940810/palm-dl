"""Aggregate experiments into manuscript tables.

Outputs (results/tables/*.csv):
  table1_cv.csv             5-fold cluster CV, all species, every test row
  table2_cv_human.csv       same folds, evaluated on the human test sites only
  table3_cv_human_only.csv  human-only training, human test
  table4_species.csv        strict cross-species holdout
  table5_hard.csv           hard-negative stratification (positives + same-protein Cys)
  table6_leakage.csv        split-protocol audit
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from config import RESULTS
from metrics import find_threshold

T = RESULTS / "tables"
T.mkdir(exist_ok=True)


def _scores(df: pd.DataFrame, thr: dict) -> tuple[np.ndarray, np.ndarray]:
    pred = np.zeros(len(df), dtype=int)
    for k, part in df.groupby("fold"):
        pred[part.index.values] = (part.p.values >= thr.get(k, 0.5)).astype(int)
    return df.y.values, pred


def pooled(df: pd.DataFrame, va: pd.DataFrame) -> dict:
    df = df.reset_index(drop=True)
    va = va.reset_index(drop=True)
    thr = {k: find_threshold(p.y.values, p.p.values) for k, p in va.groupby("fold")}
    # thresholds must be indexed by the fold values present in df
    y, pred = _scores(df, thr)
    yp, p = y, df.p.values
    tn = int(((yp == 0) & (pred == 0)).sum())
    fp = int(((yp == 0) & (pred == 1)).sum())
    fn = int(((yp == 1) & (pred == 0)).sum())
    tp = int(((yp == 1) & (pred == 1)).sum())
    denom = np.sqrt(max((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn), 1e-9))
    return {
        "n": len(df), "n_pos": int(yp.sum()),
        "auc": roc_auc_score(yp, p) if len(np.unique(yp)) > 1 else np.nan,
        "aupr": average_precision_score(yp, p),
        "mcc": (tp * tn - fp * fn) / denom,
        "acc": (tp + tn) / len(df),
        "sn": tp / max(tp + fn, 1), "sp": tn / max(tn + fp, 1),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
    }


def balanced_parts(df: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    parts = []
    for _, part in df.groupby("fold"):
        pos = part[part.y == 1]
        neg = part[part.y == 0]
        neg = neg.sample(n=min(len(pos), len(neg)), random_state=seed)
        parts.append(pd.concat([pos, neg]))
    return pd.concat(parts).reset_index(drop=True)


def table(pred_file: str, models: list[str], subset: str = "all") -> pd.DataFrame:
    path = RESULTS / pred_file
    if not path.exists():
        return pd.DataFrame()
    d = pd.read_csv(path)
    rows = []
    for m in models:
        sub = d[d.model == m]
        if sub.empty:
            continue
        te, va = sub[sub.set == "test"], sub[sub.set == "val"]
        if subset == "human" and "is_human" in sub.columns:
            te = te[te.is_human.astype(bool)]
            va = va[va.is_human.astype(bool)]
        if subset == "hard":
            te = te[(te.y == 1) | (te.source == "same_protein_cys")]
        if te.empty or va.empty:
            continue
        r = pooled(te, va)
        rb = pooled(balanced_parts(te), balanced_parts(va))
        rows.append({"model": m,
                     **{f"nat_{k}": v for k, v in r.items()},
                     **{f"bal_{k}": v for k, v in rb.items()}})
    return pd.DataFrame(rows)


def show(name: str, df: pd.DataFrame) -> None:
    if df.empty:
        print(f"\n== {name} == (missing)")
        return
    df.to_csv(T / f"{name}.csv", index=False)
    cols = [c for c in ["model", "nat_n", "nat_n_pos", "nat_auc", "nat_aupr", "nat_mcc",
                        "nat_acc", "nat_sn", "nat_sp", "bal_auc", "bal_mcc"] if c in df.columns]
    print(f"\n== {name} ==\n{df[cols].round(3).to_string(index=False)}")


def main() -> None:
    show("table1_cv", table("preds_cv_all.csv", ["pwm", "lr", "rf", "cnn", "esm", "ensemble"]))
    show("table2_cv_human", table("preds_cv_all.csv", ["pwm", "lr", "rf", "cnn", "esm", "ensemble"], "human"))
    show("table3_cv_human_only", table("preds_cvh_all.csv", ["pwm", "lr", "rf", "cnn", "esm", "ensemble"], "human"))
    # species experiment: metrics are stored per fold in metrics_species.csv (no val preds saved)
    sp = RESULTS / "metrics_species.csv"
    if sp.exists():
        d = pd.read_csv(sp)
        t4 = (d[d.protocol.isin(["test_natural", "test_balanced"])]
              .pivot_table(index="model", columns="protocol",
                           values=["mcc", "auc", "aupr", "acc", "sn", "sp", "ece"], aggfunc="mean"))
        t4.columns = [f"{a}_{b.replace('test_', '')}" for a, b in t4.columns]
        t4 = t4.reset_index()
        cols = [c for c in ["model", "mcc_natural", "auc_natural", "aupr_natural", "acc_natural",
                            "sn_natural", "sp_natural", "mcc_balanced", "auc_balanced"] if c in t4.columns]
        t4 = t4[cols].round(4)
        t4.to_csv(T / "table4_species.csv", index=False)
        print("\n== table4_species ==\n", t4.to_string(index=False))
    show("table5_hard", table("preds_cv_all.csv", ["pwm", "lr", "rf", "cnn", "esm", "ensemble"], "hard"))
    if (RESULTS / "leakage_ladder.csv").exists():
        (pd.read_csv(RESULTS / "leakage_ladder.csv")).to_csv(T / "table6_leakage.csv", index=False)
        print("\n== table6_leakage ==\n", pd.read_csv(RESULTS / "leakage_ladder.csv").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
