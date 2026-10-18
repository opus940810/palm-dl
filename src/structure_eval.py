"""Does predicted structure improve PalmDL on human sites?

Paired comparison on the subset of human sites with an AlphaFold model:
the same random forest is trained with and without three structural features
(pLDDT at the cysteine, mean pLDDT over +/-5 residues, C-alpha packing density
within 10 A), using the existing cluster-disjoint folds.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from metrics import find_threshold
from models import make_rf, set_seed

OUT = RESULTS / "structure"


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], validate="one_to_one")
    df = df[(df.window.str.len() == WINDOW) & df.is_human].reset_index(drop=True)
    st = pd.read_csv(OUT / "site_structure_features.csv")
    df = df.merge(st, on=["acc", "pos"], how="inner")
    print(f"human sites with structure: {len(df)} ({int(df.label.sum())} positives, "
          f"{df.acc.nunique()} proteins)")

    X = encode_flat(df.window.tolist())
    C = context_features(df)
    S = df[["plddt_site", "plddt_window", "packing"]].values.astype(np.float32)
    S = (S - S.mean(0)) / (S.std(0) + 1e-6)

    rows = []
    for k in sorted(df.fold.unique()):
        te = np.where(df.fold.values == k)[0]
        tr = np.where(df.fold.values != k)[0]
        if len(te) < 20 or df.label.values[te].min() == df.label.values[te].max():
            continue
        # inner validation for early-free RF is unnecessary; use fixed hyperparameters
        for name, Xtr, Xte in (("base", np.hstack([X, C]), np.hstack([X, C])),
                               ("base+structure", np.hstack([X, C, S]), np.hstack([X, C, S]))):
            set_seed(SEED + k)
            rf = make_rf()
            rf.fit(Xtr[tr], df.label.values[tr])
            p = rf.predict_proba(Xte[te])[:, 1]
            rows.append({"fold": k, "model": name, "n": len(te), "n_pos": int(df.label.values[te].sum()),
                         "auc": roc_auc_score(df.label.values[te], p)})
        print(f"fold {k}: base {rows[-2]['auc']:.3f} vs +structure {rows[-1]['auc']:.3f}", flush=True)

    res = pd.DataFrame(rows)
    piv = res.pivot(index="fold", columns="model", values="auc")
    diffs = (piv["base+structure"] - piv["base"]).values
    summary = {
        "n_sites": int(len(df)), "n_pos": int(df.label.sum()), "n_proteins": int(df.acc.nunique()),
        "base_auc": float(piv["base"].mean()), "struct_auc": float(piv["base+structure"].mean()),
        "mean_delta": float(diffs.mean()), "folds": len(piv),
        "delta_positive_folds": int((diffs > 0).sum()),
        "structural_feature_stats": {
            "plddt_site_median": float(df.plddt_site.median()),
            "plddt_pos_mean": float(df[df.label == 1].plddt_site.mean()),
            "plddt_neg_mean": float(df[df.label == 0].plddt_site.mean()),
            "packing_pos_mean": float(df[df.label == 1].packing.mean()),
            "packing_neg_mean": float(df[df.label == 0].packing.mean()),
        },
    }
    res.to_csv(OUT / "paired_cv.csv", index=False)
    (OUT / "structure_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
