"""Hybrid structure-aware PalmDL: structure where available, graceful fallback.

The previous prototype trained on the 1,974 structured sites only, so it lost the
19,000 sequence-only training sites and under-performed on the proteome scan.
Here the random forest is trained on the full fold training set with three
AlphaFold features plus a missing-value indicator; proteins without a model get
the median and flag 0, which lets the tree fall back to the sequence path.

Paired evaluation: same folds, same training rows, with and without the
structure block, on all test rows and on the structured subset.
"""
from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, MODELS, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from models import make_rf, set_seed

OUT = RESULTS / "structure"
FEATS = ["plddt_site", "plddt_window", "packing"]


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    st = pd.read_csv(OUT / "site_structure_features.csv")
    df = df.merge(st, on=["acc", "pos"], how="left")
    has = df.plddt_site.notna().values
    print(f"rows {len(df)}; with structure {int(has.sum())} "
          f"({int(df[has].label.sum())} positives, {df[has].acc.nunique()} proteins)")

    med = df.loc[has, FEATS].median()
    S = df[FEATS].fillna(med).values.astype(np.float32)
    S = np.column_stack([S, has.astype(np.float32)])
    S = (S - S.mean(0)) / (S.std(0) + 1e-6)

    X = encode_flat(df.window.tolist())
    C = context_features(df)
    base = np.hstack([X, C])
    hybrid = np.hstack([X, C, S])
    y = df.label.values

    rows = []
    for k in sorted(df.fold.unique()):
        te = np.where(df.fold.values == k)[0]
        tr = np.where(df.fold.values != k)[0]
        inner = df.iloc[tr]
        bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
        n_pos = int(inner.label.sum())
        keep = bg.sample(n=min(len(bg), n_pos), random_state=SEED + k)
        fit = np.concatenate([inner[inner.source != "background_cys"].index.values, keep.index.values])
        for name, M in (("base", base), ("hybrid", hybrid)):
            set_seed(SEED + k)
            rf = make_rf().fit(M[fit], y[fit])
            p = rf.predict_proba(M[te])[:, 1]
            sub = te[has[te]]
            rows.append({"fold": k, "model": name,
                         "auc_all": roc_auc_score(y[te], p),
                         "auc_structured": roc_auc_score(y[sub], p[has[te]]) if len(sub) > 20 else np.nan})
        print(f"fold {k}: all {rows[-2]['auc_all']:.4f}->{rows[-1]['auc_all']:.4f} | "
              f"structured {rows[-2]['auc_structured']:.4f}->{rows[-1]['auc_structured']:.4f}", flush=True)

    res = pd.DataFrame(rows)
    piv_all = res.pivot(index="fold", columns="model", values="auc_all")
    piv_str = res.pivot(index="fold", columns="model", values="auc_structured")
    summary = {
        "n_rows": int(len(df)), "n_structured": int(has.sum()),
        "auc_all_base": float(piv_all["base"].mean()),
        "auc_all_hybrid": float(piv_all["hybrid"].mean()),
        "delta_all": float((piv_all["hybrid"] - piv_all["base"]).mean()),
        "auc_struct_base": float(piv_str["base"].mean()),
        "auc_struct_hybrid": float(piv_str["hybrid"].mean()),
        "delta_struct": float((piv_str["hybrid"] - piv_str["base"]).mean()),
        "folds_positive_all": int(((piv_all["hybrid"] - piv_all["base"]) > 0).sum()),
        "folds_positive_struct": int(((piv_str["hybrid"] - piv_str["base"]) > 0).sum()),
    }
    res.to_csv(OUT / "hybrid_cv.csv", index=False)
    (OUT / "hybrid_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))

    # final artifact: train on all rows, save
    set_seed(SEED)
    rf = make_rf().fit(hybrid, y)
    joblib.dump({"model": rf, "median": med.to_dict(), "features": FEATS},
                MODELS / "rf_hybrid_struct.joblib")
    print("saved models/rf_hybrid_struct.joblib")


if __name__ == "__main__":
    main()
