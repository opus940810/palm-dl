"""Attach the training-time scaler to the hybrid model artifact.

structure_hybrid.py standardised the structure block with the mean/std of the
full training matrix, but only the model and the imputation medians were saved.
Inference must use exactly those numbers, so we recompute them from the training
table and store them alongside the model.
"""
from __future__ import annotations

import joblib
import numpy as np
import pandas as pd

from config import DATA_PROC, MODELS, RESULTS, WINDOW

FEATS = ["plddt_site", "plddt_window", "packing"]


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window"]], on=["acc", "pos"], validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    st = pd.read_csv(RESULTS / "structure" / "site_structure_features.csv")
    df = df.merge(st, on=["acc", "pos"], how="left")
    has = df.plddt_site.notna().values
    med = df.loc[has, FEATS].median()
    S = df[FEATS].fillna(med).values.astype(np.float32)
    S = np.column_stack([S, has.astype(np.float32)])
    mean, std = S.mean(0), S.std(0) + 1e-6
    path = MODELS / "rf_hybrid_struct.joblib"
    obj = joblib.load(path)
    obj["scaler_mean"] = mean
    obj["scaler_std"] = std
    obj["median"] = med.to_dict()
    obj["features"] = FEATS
    joblib.dump(obj, path)
    print("scaler mean", np.round(mean, 3))
    print("scaler std ", np.round(std, 3))
    print("saved", path)


if __name__ == "__main__":
    main()
