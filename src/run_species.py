"""Independent cross-species evaluation.

Training uses non-human sites only; every annotation that was propagated from a
human entry has already been removed by `splits.py`, and any protein cluster
containing a human test protein is excluded from training. Human background
negatives (proteins without any lipidation annotation) are split 50/50 by
protein hash so that the test negatives are protein-disjoint from training.
"""
from __future__ import annotations

import argparse
import hashlib
import time

import numpy as np
import pandas as pd
import torch

from config import DATA_PROC, RESULTS, SEED, WINDOW
from features import context_features, encode_flat, encode_onehot
from metrics import find_threshold, metrics_at
from models import PalmCNN, make_logreg, make_rf, predict_torch, set_seed, train_torch
from run_cv import balanced_subset, pwm_scores

torch.set_num_threads(12)


def half_of(acc: str) -> int:
    return int(hashlib.md5(acc.encode()).hexdigest(), 16) % 2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="pwm,lr,cnn")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--bg-ratio", type=float, default=1.0)
    ap.add_argument("--tag", default="species")
    args = ap.parse_args()

    set_seed(SEED)
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    sp = pd.read_csv(DATA_PROC / "split_species.csv")
    df = sp.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                  on=["acc", "pos"], how="left", validate="one_to_one")
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    df["half"] = df.acc.map(half_of)

    is_test = (df.split_species == "test")
    # human background negatives: half for training, half kept for testing
    bg_human = is_test & (df.source == "background_cys")
    train_mask = (df.split_species == "train") | (bg_human & (df.half == 0))
    test_mask = is_test & ~(bg_human & (df.half == 0))
    train_idx = np.where(train_mask.values)[0]
    test_idx = np.where(test_mask.values)[0]
    print(f"train rows={len(train_idx)} (pos={int(df.iloc[train_idx].label.sum())}), "
          f"test rows={len(test_idx)} (pos={int(df.iloc[test_idx].label.sum())})", flush=True)

    X1h = encode_onehot(df.window.tolist())
    Cctx = context_features(df)
    Xflat = encode_flat(df.window.tolist())

    tr = df.iloc[train_idx]
    clusters = tr.cluster.unique()
    rng = np.random.default_rng(SEED)
    val_clusters = set(rng.choice(clusters, size=max(2, int(0.15 * len(clusters))), replace=False))
    is_val = tr.cluster.isin(val_clusters).values
    val_idx = train_idx[is_val]
    inner_idx = train_idx[~is_val]
    inner = df.iloc[inner_idx]
    bg = inner[(inner.label == 0) & (inner.source == "background_cys")]
    n_pos = int(inner.label.sum())
    keep_bg = bg.sample(n=min(len(bg), int(args.bg_ratio * n_pos)), random_state=SEED)
    fit_idx = np.concatenate([inner[inner.source != "background_cys"].index.values, keep_bg.index.values])
    print(f"fit rows={len(fit_idx)} (pos={int(df.iloc[fit_idx].label.sum())}), val rows={len(val_idx)}", flush=True)

    y_val = df.label.values[val_idx]
    bal_val = balanced_subset(df.iloc[val_idx], np.zeros(len(val_idx)), SEED)
    bal_test = balanced_subset(df.iloc[test_idx], np.zeros(len(test_idx)), SEED + 7)
    rows, preds = [], []

    for name in [m.strip() for m in args.models.split(",")]:
        t0 = time.time()
        set_seed(SEED)
        if name == "pwm":
            p_val = pwm_scores(df.iloc[fit_idx], df.iloc[val_idx])
            p_test = pwm_scores(df.iloc[fit_idx], df.iloc[test_idx])
        elif name == "lr":
            clf = make_logreg()
            clf.fit(Xflat[fit_idx], df.label.values[fit_idx])
            p_val = clf.predict_proba(Xflat[val_idx])[:, 1]
            p_test = clf.predict_proba(Xflat[test_idx])[:, 1]
        elif name == "rf":
            clf = make_rf()
            clf.fit(np.hstack([Xflat[fit_idx], Cctx[fit_idx]]), df.label.values[fit_idx])
            p_val = clf.predict_proba(np.hstack([Xflat[val_idx], Cctx[val_idx]]))[:, 1]
            p_test = clf.predict_proba(np.hstack([Xflat[test_idx], Cctx[test_idx]]))[:, 1]
        elif name == "cnn":
            net = PalmCNN()
            net, best, hist = train_torch(
                net, X1h[fit_idx], Cctx[fit_idx], df.label.values[fit_idx].astype(np.float32),
                X1h[val_idx], Cctx[val_idx], y_val.astype(np.float32),
                epochs=args.epochs)
            print(f"  cnn stopped at epoch {best['epoch']} (val_auc={best['score']:.4f})", flush=True)
            p_val = predict_torch(net, X1h[val_idx], Cctx[val_idx])
            p_test = predict_torch(net, X1h[test_idx], Cctx[test_idx])
        else:
            raise ValueError(name)

        t_nat = find_threshold(y_val, p_val)
        t_bal = find_threshold(y_val[bal_val], p_val[bal_val])
        m_nat = metrics_at(df.label.values[test_idx], p_test, t_nat)
        m_bal = metrics_at(df.label.values[test_idx][bal_test], p_test[bal_test], t_bal)
        m_val = metrics_at(y_val, p_val, t_nat)
        for tag, m in (("test_natural", m_nat), ("test_balanced", m_bal), ("val_natural", m_val)):
            rows.append({"model": name, "protocol": tag, "t_nat": t_nat, "t_bal": t_bal,
                         **m, "sec": round(time.time() - t0, 1)})
        preds.append(pd.DataFrame({"model": name, "acc": df.acc.values[test_idx], "pos": df.pos.values[test_idx],
                                   "y": df.label.values[test_idx], "p": p_test,
                                   "source": df.source.values[test_idx]}))
        print(f"  {name:4s} human test: MCC={m_nat['mcc']:.3f} AUC={m_nat['auc']:.3f} "
              f"AUPR={m_nat['aupr']:.3f} | balanced MCC={m_bal['mcc']:.3f}", flush=True)

    pd.concat(preds).to_csv(RESULTS / f"preds_{args.tag}.csv", index=False)
    pd.DataFrame(rows).to_csv(RESULTS / f"metrics_{args.tag}.csv", index=False)
    print(pd.DataFrame(rows)[["model", "protocol", "mcc", "auc", "aupr", "acc", "sn", "sp", "ece"]].round(4).to_string())


if __name__ == "__main__":
    main()
