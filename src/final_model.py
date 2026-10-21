"""Train the deployable PalmDL model on all labelled sites and score new sites.

The cross-validated ensemble in the paper averages within-fold ranks.  A released
model cannot rank against a fold, so the deployable version averages logits of
the two calibrated views and then applies Platt scaling; on the same folds this
changes AUC by <0.003 (see results/headline_stats.json).
"""
from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

import esm
from config import AA, DATA_PROC, MODELS, ROOT, SEED, WINDOW
from features import context_features, encode_flat, encode_onehot
from metrics import find_threshold
from models import EmbeddingMLP, make_rf, set_seed, train_torch
from run_cv import balanced_subset

torch.set_num_threads(4)
CKPT = ROOT.parent / ".cache" / "esm2_t12_35M_UR50D.pt"
RF_PATH = MODELS / "rf_final.joblib"
ESM_PATH = MODELS / "esm_final.pt"
CAL_PATH = MODELS / "calibrator.joblib"
META_PATH = MODELS / "final_meta.json"


def load_esm(path=CKPT):
    orig = torch.load

    def _load(*a, **k):
        k["weights_only"] = False
        return orig(*a, **k)

    torch.load = _load
    try:
        return esm.pretrained.load_model_and_alphabet_local(str(path))
    finally:
        torch.load = orig


_ESM_CACHE = {}


def embed_windows(windows: list[str], bs: int = 64) -> np.ndarray:
    """Mean-pooled + centre-residue ESM-2 representations for arbitrary windows."""
    if "model" not in _ESM_CACHE:
        model, alphabet = load_esm()
        model.eval()
        _ESM_CACHE.update(model=model, bc=alphabet.get_batch_converter(), dim=model.embed_dim)
    model, bc, dim = _ESM_CACHE["model"], _ESM_CACHE["bc"], _ESM_CACHE["dim"]
    import re
    windows = [re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", w) for w in windows]
    out = np.zeros((len(windows), 2 * dim), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(windows), bs):
            chunk = windows[start:start + bs]
            _, _, toks = bc([(f"s{i}", w) for i, w in enumerate(chunk)])
            rep = model(toks, repr_layers=[model.num_layers], return_contacts=False)[
                "representations"][model.num_layers]
            out[start:start + len(chunk), :dim] = rep[:, 1:-1].mean(1).numpy()
            out[start:start + len(chunk), dim:] = rep[:, 1 + WINDOW // 2].numpy()
    return out


def logit(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))


def train() -> dict:
    set_seed(SEED)
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    # 66 UniProt features map to a non-cysteine residue; exclude them from the
    # deployed model (their removal does not change CV performance, Table S1)
    keep = ~((df.label == 1) & (df.window.str[15] != "C")).values
    df = df[keep].reset_index(drop=True)
    Xflat = encode_flat(df.window.tolist())
    Cctx = context_features(df)
    Xesm = np.hstack([np.load(DATA_PROC.parent / "embeddings" / "esm2_mean.npy").astype(np.float32),
                      np.load(DATA_PROC.parent / "embeddings" / "esm2_center.npy").astype(np.float32)])[keep]

    rng = np.random.default_rng(SEED)
    val_clusters = set(rng.choice(df.cluster.unique(), size=max(2, int(0.15 * df.cluster.nunique())),
                                  replace=False))
    is_val = df.cluster.isin(val_clusters).values
    val_idx = np.where(is_val)[0]
    inner = df[~is_val]
    bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
    n_pos = int(inner.label.sum())
    keep = bg.sample(n=min(len(bg), n_pos), random_state=SEED)
    fit_idx = np.concatenate([inner[inner.source != "background_cys"].index.values, keep.index.values])
    y = df.label.values.astype(np.float32)
    print(f"fit={len(fit_idx)} val={len(val_idx)} positives={int(y.sum())}")

    rf = make_rf()
    rf.fit(np.hstack([Xflat[fit_idx], Cctx[fit_idx]]), y[fit_idx])
    joblib.dump(rf, RF_PATH)

    net = EmbeddingMLP(in_dim=Xesm.shape[1])
    net, best, _ = train_torch(net, Xesm[fit_idx], Cctx[fit_idx], y[fit_idx],
                               Xesm[val_idx], Cctx[val_idx], y[val_idx], epochs=40, lr=5e-4)
    torch.save(net.state_dict(), ESM_PATH)

    def prob(model, X, C):
        model.eval()
        with torch.no_grad():
            return torch.sigmoid(model(torch.as_tensor(X), torch.as_tensor(C))).numpy()

    p_rf_v = rf.predict_proba(np.hstack([Xflat[val_idx], Cctx[val_idx]]))[:, 1]
    p_esm_v = prob(net, Xesm[val_idx], Cctx[val_idx])
    z_val = 0.5 * (logit(p_rf_v) + logit(p_esm_v))
    cal = LogisticRegression(C=1.0, max_iter=1000).fit(z_val.reshape(-1, 1), y[val_idx])
    joblib.dump(cal, CAL_PATH)
    p_cal = cal.predict_proba(z_val.reshape(-1, 1))[:, 1]

    bal = balanced_subset(df.iloc[val_idx], p_cal, SEED)
    meta = {
        "val_auc": float(roc_auc_score(y[val_idx], p_cal)),
        "val_auc_rf": float(roc_auc_score(y[val_idx], p_rf_v)),
        "val_auc_esm": float(roc_auc_score(y[val_idx], p_esm_v)),
        "threshold_natural": find_threshold(y[val_idx], p_cal),
        "threshold_balanced": find_threshold(y[val_idx][bal], p_cal[bal]),
        "n_fit": int(len(fit_idx)), "n_val": int(len(val_idx)),
        "window": WINDOW, "esm": "esm2_t12_35M_UR50D", "rf_trees": 400,
    }
    META_PATH.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))
    return meta


def load():
    rf = joblib.load(RF_PATH)
    net = EmbeddingMLP(in_dim=960)
    net.load_state_dict(torch.load(ESM_PATH, map_location="cpu"))
    net.eval()
    cal = joblib.load(CAL_PATH)
    meta = json.loads(META_PATH.read_text(encoding="utf-8"))
    return rf, net, cal, meta


def score(df: pd.DataFrame, batch: bool = True) -> pd.DataFrame:
    """Score a dataframe with columns acc, pos, window, compartment, protein_length."""
    rf, net, cal, meta = load()
    df = df.copy().reset_index(drop=True)
    Xflat = encode_flat(df.window.tolist())
    Cctx = context_features(df)
    Xesm = embed_windows(df.window.tolist())
    p_rf = rf.predict_proba(np.hstack([Xflat, Cctx]))[:, 1]
    with torch.no_grad():
        p_esm = torch.sigmoid(net(torch.as_tensor(Xesm), torch.as_tensor(Cctx))).numpy()
    z = 0.5 * (logit(p_rf) + logit(p_esm))
    p = cal.predict_proba(z.reshape(-1, 1))[:, 1]
    return df.assign(p_rf=p_rf, p_esm=p_esm, p=p)


if __name__ == "__main__":
    train()
