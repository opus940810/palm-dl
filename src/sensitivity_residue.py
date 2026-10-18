"""Sensitivity analysis: drop annotations whose canonical residue is not C."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from metrics import find_threshold, metrics_at
from models import make_rf, set_seed
from run_cv import balanced_subset


def run(df: pd.DataFrame, label: str) -> dict:
    X = encode_flat(df.window.tolist())
    ctx = context_features(df)
    aucs, mccs = [], []
    for k in sorted(df.fold.unique()):
        te = np.where(df.fold.values == k)[0]
        tr = np.where(df.fold.values != k)[0]
        tr_clusters = df.iloc[tr].cluster.unique()
        rng = np.random.default_rng(SEED + k)
        vc = set(rng.choice(tr_clusters, size=max(2, int(0.15 * len(tr_clusters))), replace=False))
        is_val = df.iloc[tr].cluster.isin(vc).values
        val = tr[is_val]
        inner = df.iloc[tr[~is_val]]
        bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
        n_pos = int(inner.label.sum())
        keep = bg.sample(n=min(len(bg), n_pos), random_state=SEED + k)
        fit = np.concatenate([inner[inner.source != "background_cys"].index.values, keep.index.values])
        set_seed(SEED + k)
        clf = make_rf()
        clf.fit(np.hstack([X[fit], ctx[fit]]), df.label.values[fit])
        p_val = clf.predict_proba(np.hstack([X[val], ctx[val]]))[:, 1]
        p_te = clf.predict_proba(np.hstack([X[te], ctx[te]]))[:, 1]
        t = find_threshold(df.label.values[val], p_val)
        b = balanced_subset(df.iloc[te], p_te, SEED + k)
        tb = find_threshold(df.label.values[val][balanced_subset(df.iloc[val], p_val, SEED + k)],
                            p_val[balanced_subset(df.iloc[val], p_val, SEED + k)])
        aucs.append(roc_auc_score(df.label.values[te], p_te))
        mccs.append(metrics_at(df.label.values[te][b], p_te[b], tb)["mcc"])
    print(f"{label}: AUC {np.mean(aucs):.4f} +/- {np.std(aucs):.4f}, balanced MCC {np.mean(mccs):.4f}")
    return {"label": label, "auc": float(np.mean(aucs)), "auc_sd": float(np.std(aucs)),
            "bal_mcc": float(np.mean(mccs)), "n": len(df), "n_pos": int(df.label.sum())}


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    splits = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = splits.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                      on=["acc", "pos"], validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    out = [run(df, "all annotations")]
    clean = df[~((df.label == 1) & (df.window.str[15] != "C"))].reset_index(drop=True)
    out.append(run(clean, "non-cysteine annotations removed"))
    pd.DataFrame(out).to_csv(RESULTS / "sensitivity_residue_filter.csv", index=False)


if __name__ == "__main__":
    main()
