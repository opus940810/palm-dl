"""Honest within-species diagonal for the cross-species transfer matrix.

The full matrix (`species_matrix.csv`) trains a model on all rows of a source
species, so its diagonal entry is optimistic.  Here each species' clusters are
split 80/20 and the diagonal is recomputed on the held-out 20 % only, using the
same RF + ESM ensemble as the off-diagonal cells.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, EMB, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from final_model import logit
from models import EmbeddingMLP, make_rf, set_seed, train_torch

torch.set_num_threads(8)
ORG = {"9606": "Human", "10090": "Mouse", "10116": "Rat", "9913": "Cattle",
       "83333": "E. coli", "3702": "A. thaliana", "9823": "Pig", "224308": "B. subtilis"}


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], validate="one_to_one")
    df["organism"] = df["organism"].astype(str).str.split(".").str[0]
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    Xflat = encode_flat(df.window.tolist())
    Cctx = context_features(df)
    Xesm = np.hstack([np.load(EMB / "esm2_mean.npy").astype(np.float32),
                      np.load(EMB / "esm2_center.npy").astype(np.float32)])

    rows = []
    for sp, name in ORG.items():
        idx = np.where((df.organism == sp).values)[0]
        if len(idx) < 40 or df.label.values[idx].sum() < 20:
            continue
        clusters = df.iloc[idx].cluster.unique()
        rng = np.random.default_rng(SEED)
        hold = set(rng.choice(clusters, size=max(2, int(0.2 * len(clusters))), replace=False))
        te = idx[df.iloc[idx].cluster.isin(hold).values]
        rest = idx[~df.iloc[idx].cluster.isin(hold).values]

        pos = df.iloc[rest][df.iloc[rest].label == 1]
        negs = df.iloc[rest][df.iloc[rest].label == 0]
        bg = df[(df.source == "background_cys") & ~df.cluster.isin(hold)]
        keep_bg = bg.sample(n=min(len(bg), len(pos)), random_state=SEED)
        fit = np.concatenate([pos.index.values, negs.index.values, keep_bg.index.values])

        left = df.iloc[rest].cluster.unique()
        val_c = set(rng.choice(left, size=max(2, int(0.15 * len(left))), replace=False))
        val = rest[df.iloc[rest].cluster.isin(val_c).values]
        fit = np.setdiff1d(fit, val)

        set_seed(SEED)
        rf = make_rf()
        rf.fit(np.hstack([Xflat[fit], Cctx[fit]]), df.label.values[fit])
        net = EmbeddingMLP(in_dim=Xesm.shape[1])
        net, _, _ = train_torch(net, Xesm[fit], Cctx[fit], df.label.values[fit].astype(np.float32),
                               Xesm[val], Cctx[val], df.label.values[val].astype(np.float32),
                               epochs=30, lr=5e-4, patience=6)
        y = df.label.values[te]
        p_rf = rf.predict_proba(np.hstack([Xflat[te], Cctx[te]]))[:, 1]
        net.eval()
        with torch.no_grad():
            p_esm = torch.sigmoid(net(torch.as_tensor(Xesm[te]),
                                      torch.as_tensor(Cctx[te]))).numpy()
        p = 1 / (1 + np.exp(-0.5 * (logit(p_rf) + logit(p_esm))))
        rows.append({"species": name, "n_test": len(te), "n_pos": int(y.sum()),
                     "auc": roc_auc_score(y, p)})
        print(rows[-1], flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "species_matrix_diagonal.csv", index=False)
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
