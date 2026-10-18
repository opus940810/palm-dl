"""Cross-validation runner for all models (PWM, logistic regression, random
forest, PalmDL-CNN).

For every outer fold the training clusters are split once more into an inner
train/validation set (cluster-disjoint). The inner validation set is used for
early stopping and for choosing the decision threshold, so the outer test fold
is never touched during model selection.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from config import DATA_PROC, EMB, RESULTS, SEED, WINDOW
from features import context_features, encode_flat, encode_onehot
from metrics import find_threshold, metrics_at
from models import EmbeddingMLP, PalmCNN, make_logreg, make_rf, predict_torch, set_seed, train_torch

import os
torch.set_num_threads(int(os.environ.get("PALMDL_THREADS", "4")))


def pwm_scores(train: pd.DataFrame, test: pd.DataFrame, eps: float = 0.5) -> np.ndarray:
    """CSS-Palm-style position-specific log-odds score of the 31-aa window."""
    from config import AA

    idx = {a: i for i, a in enumerate(AA)}
    fp = np.full((WINDOW, 20), eps)
    fn = np.full((WINDOW, 20), eps)
    for w, y in zip(train.window, train.label):
        mat = fp if y == 1 else fn
        for j, ch in enumerate(w):
            if ch in idx:
                mat[j, idx[ch]] += 1
    fp /= fp.sum(axis=1, keepdims=True)
    fn /= fn.sum(axis=1, keepdims=True)
    lo = np.log2(fp / fn)
    S = np.zeros(len(test))
    for i, w in enumerate(test.window):
        s = 0.0
        for j, ch in enumerate(w):
            if ch in idx:
                s += lo[j, idx[ch]]
        S[i] = s
    # map to (0,1) with a logistic so that thresholds are comparable
    return 1.0 / (1.0 + np.exp(-S / 5.0))


def balanced_subset(df: pd.DataFrame, p: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    pos = np.where(df.label.values == 1)[0]
    neg = np.where(df.label.values == 0)[0]
    neg = rng.choice(neg, size=min(len(pos), len(neg)), replace=False)
    return np.concatenate([pos, neg])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="pwm,lr,rf,cnn")
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--bg-ratio", type=float, default=1.0,
                    help="background negatives kept per positive in training")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--subset", default="all", choices=["all", "human"],
                    help="restrict the whole CV to human proteins")
    ap.add_argument("--tag", default="cv")
    args = ap.parse_args()

    set_seed(SEED)
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], how="left", validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    human_mask = df.is_human.astype(bool).values
    if args.subset == "human":
        df = df[human_mask].reset_index(drop=True)
    print(f"dataset: {len(df)} rows, {int(df.label.sum())} positives")

    X1h = encode_onehot(df.window.tolist())
    Cctx = context_features(df)
    Xflat = encode_flat(df.window.tolist())
    print(f"features: onehot {X1h.shape}, flat {Xflat.shape}, ctx {Cctx.shape}")

    models = [m.strip() for m in args.models.split(",")]
    Xesm = None
    if "esm" in models:
        Xesm = np.hstack([np.load(EMB / "esm2_mean.npy").astype(np.float32),
                          np.load(EMB / "esm2_center.npy").astype(np.float32)])
        if len(Xesm) == len(human_mask):
            Xesm = Xesm[human_mask] if args.subset == "human" else Xesm
        if len(Xesm) != len(df):
            raise SystemExit(f"embedding rows {len(Xesm)} != dataset rows {len(df)}; re-run esm_embed.py")
        print(f"esm features: {Xesm.shape}")
    folds = [int(f) for f in args.folds.split(",")]
    preds, rows_hist = [], []

    for k in folds:
        t0 = time.time()
        test_idx = np.where(df.fold.values == k)[0]
        train_idx = np.where(df.fold.values != k)[0]
        tr_clusters = df.iloc[train_idx].cluster.unique()
        rng = np.random.default_rng(SEED + k)
        val_clusters = set(rng.choice(tr_clusters, size=max(2, int(0.15 * len(tr_clusters))), replace=False))
        is_val = df.iloc[train_idx].cluster.isin(val_clusters).values
        val_idx = train_idx[is_val]
        inner_idx = train_idx[~is_val]

        inner = df.iloc[inner_idx]
        bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
        n_pos = int(inner.label.sum())
        keep_bg = bg.sample(n=min(len(bg), int(args.bg_ratio * n_pos)), random_state=SEED + k)
        fit_idx = np.concatenate([inner[inner.source != "background_cys"].index.values, keep_bg.index.values])
        print(f"[fold {k}] fit={len(fit_idx)} (pos={int(df.iloc[fit_idx].label.sum())}) "
              f"val={len(val_idx)} test={len(test_idx)}", flush=True)

        y_val = df.label.values[val_idx]
        bal_val = balanced_subset(df.iloc[val_idx], np.zeros(len(val_idx)), SEED + k)
        yv_bal = y_val[bal_val]

        for name in models:
            set_seed(SEED + k)
            t1 = time.time()
            if name == "pwm":
                p_val = pwm_scores(df.iloc[fit_idx], df.iloc[val_idx])
                p_test = pwm_scores(df.iloc[fit_idx], df.iloc[test_idx])
                hist = []
            elif name == "lr":
                clf = make_logreg()
                clf.fit(Xflat[fit_idx], df.label.values[fit_idx])
                p_val = clf.predict_proba(Xflat[val_idx])[:, 1]
                p_test = clf.predict_proba(Xflat[test_idx])[:, 1]
                hist = []
            elif name == "rf":
                clf = make_rf()
                clf.fit(np.hstack([Xflat[fit_idx], Cctx[fit_idx]]), df.label.values[fit_idx])
                p_val = clf.predict_proba(np.hstack([Xflat[val_idx], Cctx[val_idx]]))[:, 1]
                p_test = clf.predict_proba(np.hstack([Xflat[test_idx], Cctx[test_idx]]))[:, 1]
                hist = []
            elif name == "cnn":
                net = PalmCNN()
                net, best, hist = train_torch(
                    net, X1h[fit_idx], Cctx[fit_idx], df.label.values[fit_idx].astype(np.float32),
                    X1h[val_idx], Cctx[val_idx], y_val.astype(np.float32),
                    epochs=args.epochs, verbose=False)
                print(f"    cnn stopped at epoch {best['epoch']} (val_auc={best['score']:.4f})", flush=True)
                p_val = predict_torch(net, X1h[val_idx], Cctx[val_idx])
                p_test = predict_torch(net, X1h[test_idx], Cctx[test_idx])
            elif name == "esm":
                net = EmbeddingMLP(in_dim=Xesm.shape[1])
                net, best, hist = train_torch(
                    net, Xesm[fit_idx], Cctx[fit_idx], df.label.values[fit_idx].astype(np.float32),
                    Xesm[val_idx], Cctx[val_idx], y_val.astype(np.float32),
                    epochs=args.epochs, lr=5e-4, verbose=False)
                print(f"    esm stopped at epoch {best['epoch']} (val_auc={best['score']:.4f})", flush=True)
                p_val = predict_torch(net, Xesm[val_idx], Cctx[val_idx])
                p_test = predict_torch(net, Xesm[test_idx], Cctx[test_idx])
            else:
                raise ValueError(name)

            t_nat = find_threshold(y_val, p_val)
            t_bal = find_threshold(yv_bal, p_val[bal_val])
            for split_name, idx, p in (("val", val_idx, p_val), ("test", test_idx, p_test)):
                preds.append(pd.DataFrame({
                    "model": name, "fold": k, "set": split_name,
                    "acc": df.acc.values[idx], "pos": df.pos.values[idx],
                    "y": df.label.values[idx], "p": p, "source": df.source.values[idx],
                    "is_human": df.is_human.values[idx],
                }))
            m_nat = metrics_at(df.label.values[test_idx], p_test, t_nat)
            bi = balanced_subset(df.iloc[test_idx], p_test, SEED + 100 + k)
            m_bal = metrics_at(df.label.values[test_idx][bi], p_test[bi], t_bal)
            m_val = metrics_at(y_val, p_val, t_nat)
            for tag, m in (("test_natural", m_nat), ("test_balanced", m_bal), ("val_natural", m_val)):
                rows_hist.append({"model": name, "fold": k, "protocol": tag,
                                  "t_nat": t_nat, "t_bal": t_bal, **m,
                                  "sec": round(time.time() - t1, 1)})
            print(f"    {name:4s} test_nat MCC={m_nat['mcc']:.3f} AUC={m_nat['auc']:.3f} "
                  f"AUPR={m_nat['aupr']:.3f} | balanced MCC={m_bal['mcc']:.3f}", flush=True)

    pd.concat(preds).to_csv(RESULTS / f"preds_{args.tag}.csv", index=False)
    met = pd.DataFrame(rows_hist)
    met.to_csv(RESULTS / f"metrics_{args.tag}.csv", index=False)
    summary = (met[met.protocol.isin(["test_natural", "test_balanced"])]
               .groupby(["model", "protocol"])[["mcc", "auc", "aupr", "acc", "sn", "sp", "f1", "ece"]]
               .agg(["mean", "std"]).round(4))
    summary.to_csv(RESULTS / f"summary_{args.tag}.csv")
    print("\n", summary.to_string())


if __name__ == "__main__":
    main()
