"""Feature encoders shared by all models."""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import AA, AA_IDX, COMPARTMENTS, FLEX, HALF, HYDROPATHY, VOLUME, WINDOW


def encode_onehot(windows: list[str], physchem: bool = True) -> np.ndarray:
    """(N, C, L) float32; channel order = 20 AA + [hydropathy, volume, flexibility]."""
    n = len(windows)
    L = WINDOW
    C = 20 + (3 if physchem else 0)
    X = np.zeros((n, C, L), dtype=np.float32)
    vol = np.array([VOLUME[a] for a in AA], dtype=np.float32)
    vol = (vol - vol.mean()) / (vol.std() + 1e-6)
    for i, w in enumerate(windows):
        assert len(w) == L, (len(w), w)
        for j, ch in enumerate(w):
            k = AA_IDX.get(ch)
            if k is not None:
                X[i, k, j] = 1.0
            if physchem:
                X[i, 20, j] = HYDROPATHY.get(ch, 0.0) / 4.5
                X[i, 21, j] = vol[AA_IDX.get(ch, 0)] if k is not None else 0.0
                X[i, 22, j] = FLEX.get(ch, 0.44)
    return X


def encode_flat(windows: list[str]) -> np.ndarray:
    """Flattened one-hot, for linear / tree baselines."""
    return encode_onehot(windows, physchem=False).reshape(len(windows), -1)


def context_features(df: pd.DataFrame) -> np.ndarray:
    """Protein-level context available at prediction time (no label leakage).

    Columns: 6 compartment one-hot, has_annotation flag, log length,
             relative position, relative distance to the closest terminus.
    """
    n = len(df)
    F = np.zeros((n, 11), dtype=np.float32)
    comp = df["compartment"].fillna("unknown").astype(str).str.lower()
    for i, c in enumerate(COMPARTMENTS):
        F[:, i] = (comp == c).astype(np.float32)
    F[:, 6] = (comp != "unknown").astype(np.float32)
    L = df["protein_length"].astype(float).clip(lower=1).values
    pos = df["pos"].astype(float).values
    F[:, 7] = np.log1p(L) / 10.0
    F[:, 8] = pos / L
    F[:, 9] = np.minimum(pos, L - pos + 1) / L
    F[:, 10] = (pos <= HALF + 1).astype(np.float32)
    return F


def load_dataset(path, split_col: str | None = None) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"organism": str})
    df["window"] = df["window"].fillna("")
    return df
