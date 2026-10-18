"""Fit a Platt calibrator for the structure-aware path.

OOF protocol: within each fold, train the hybrid RF on the fold's training rows
(positives + same-protein negatives + 1:1 background), predict the held-out
fold, combine with the cached ESM out-of-fold view by logit averaging, then fit
one logistic calibrator on all out-of-fold combinations.
"""
from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, MODELS, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from final_model import logit
from models import make_rf, set_seed

FEATS = ["plddt_site", "plddt_window", "packing"]


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    st = pd.read_csv(RESULTS / "structure" / "site_structure_features.csv")
    df = df.merge(st, on=["acc", "pos"], how="left")
    obj = joblib.load(MODELS / "rf_hybrid_struct.joblib")
    med, mean, std = obj["median"], obj["scaler_mean"], obj["scaler_std"]

    has = df.plddt_site.notna().values
    S = df[FEATS].fillna(pd.Series(med)).values.astype(np.float32)
    S = np.column_stack([S, has.astype(np.float32)])
    S = (S - mean) / std
    X = encode_flat(df.window.tolist())
    C = context_features(df)
    H = np.hstack([X, C, S])
    y = df.label.values

    esm = pd.read_csv(RESULTS / "preds_cv_esm.csv")
    esm = esm[(esm.model == "esm") & (esm.set == "test")][["acc", "pos", "fold", "p"]]
    esm = esm.rename(columns={"p": "p_esm"})

    parts = []
    for k in sorted(df.fold.unique()):
        te = np.where(df.fold.values == k)[0]
        tr = np.where(df.fold.values != k)[0]
        inner = df.iloc[tr]
        bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
        n_pos = int(inner.label.sum())
        keep = bg.sample(n=min(len(bg), n_pos), random_state=SEED + k)
        fit = np.concatenate([inner[inner.source != "background_cys"].index.values, keep.index.values])
        set_seed(SEED + k)
        rf = make_rf().fit(H[fit], y[fit])
        p = rf.predict_proba(H[te])[:, 1]
        part = df.iloc[te][["acc", "pos", "fold", "label"]].copy()
        part["p_hyb"] = p
        parts.append(part)
        print(f"fold {k}: hybrid RF AUC {roc_auc_score(y[te], p):.4f}", flush=True)

    oof = pd.concat(parts, ignore_index=True)
    oof = oof.merge(esm, on=["acc", "pos", "fold"], how="inner")
    z = 0.5 * (logit(oof.p_hyb.values) + logit(oof.p_esm.values))
    cal = LogisticRegression(C=1.0, max_iter=1000).fit(z.reshape(-1, 1), oof.label.values)
    p_cal = cal.predict_proba(z.reshape(-1, 1))[:, 1]
    meta = {"n_oof": int(len(oof)),
            "auc_hybrid_rf": float(roc_auc_score(oof.label, oof.p_hyb)),
            "auc_hybrid_ensemble": float(roc_auc_score(oof.label, p_cal))}
    joblib.dump(cal, MODELS / "calibrator_hybrid.joblib")
    (MODELS / "calibrator_hybrid_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    oof.assign(z=z, p_cal=p_cal).to_csv(RESULTS / "structure" / "hybrid_oof.csv", index=False)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
