"""Apply the hybrid structure-aware RF to the scan subset and the literature set."""
from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, MODELS, RESULTS
from features import context_features, encode_flat

OUT = RESULTS / "structure"
FEATS = ["plddt_site", "plddt_window", "packing"]
CACHE = MODELS / "rf_hybrid_struct.joblib"


def predict(frame: pd.DataFrame) -> pd.DataFrame:
    rf_base = joblib.load(MODELS / "rf_final.joblib")
    obj = joblib.load(CACHE)
    rf_hyb, med = obj["model"], obj["median"]
    mean, std = obj["scaler_mean"], obj["scaler_std"]
    X = encode_flat(frame.window.tolist())
    C = context_features(frame)
    S = frame[FEATS].astype(float).fillna(pd.Series(med))
    has = frame.plddt_site.notna().astype(float).values.reshape(-1, 1)
    S = np.column_stack([S.values, has])
    S = (S - mean) / std
    out = frame.copy()
    out["p_base"] = rf_base.predict_proba(np.hstack([X, C]))[:, 1]
    out["p_hybrid"] = rf_hyb.predict_proba(np.hstack([X, C, S]))[:, 1]
    return out


def main() -> None:
    res = {}

    # scan subset
    scan = pd.read_csv(OUT / "scan_struct_scored.csv")
    d = predict(scan)
    known = dict(zip(zip(d.acc, d.pos), d.known_uniprot.fillna(False)))
    for k in (322, 668):
        a = set(d.nlargest(k, "p_base")[["acc", "pos"]].itertuples(index=False, name=None))
        b = set(d.nlargest(k, "p_hybrid")[["acc", "pos"]].itertuples(index=False, name=None))
        res[f"scan_top{k}_known_base"] = int(sum(known[x] for x in a))
        res[f"scan_top{k}_known_hybrid"] = int(sum(known[x] for x in b))
        res[f"scan_top{k}_overlap"] = int(len(a & b))
        print(f"scan top {k}: base {res[f'scan_top{k}_known_base']} -> hybrid "
              f"{res[f'scan_top{k}_known_hybrid']}, overlap {res[f'scan_top{k}_overlap']}")

    # literature set
    lit = pd.read_csv(OUT / "literature_struct_scored.csv")
    l = predict(lit)
    if l.label.nunique() == 2:
        res["lit_auc_base_rf"] = float(roc_auc_score(l.label, l.p_base))
        res["lit_auc_hybrid_rf"] = float(roc_auc_score(l.label, l.p_hybrid))
        print(f"literature RF: base {res['lit_auc_base_rf']:.3f} -> hybrid "
              f"{res['lit_auc_hybrid_rf']:.3f}")
    res["n_scan"] = int(len(d))
    res["n_lit"] = int(len(l))
    (OUT / "hybrid_apply.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
