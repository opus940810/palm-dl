"""Cross-species transfer matrix.

For each source species with enough annotated sites we train PalmDL (RF + ESM
views) on that species only and evaluate on every target species.  The diagonal
is therefore within-species performance, off-diagonal entries measure how much
of the palmitoylation code transfers between organisms.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, EMB, MODELS, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from final_model import logit
from metrics import find_threshold, metrics_at
from models import EmbeddingMLP, make_rf, set_seed, train_torch

torch.set_num_threads(4)
ORG = {"9606": "Human", "10090": "Mouse", "10116": "Rat", "9913": "Cattle",
       "83333": "E. coli", "3702": "A. thaliana", "9823": "Pig", "224308": "B. subtilis"}


def main() -> None:
    set_seed(SEED)
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

    counts = df[df.label == 1].organism.value_counts()
    species = [s for s in ORG if s in counts.index and counts[s] >= 75]
    print("species:", [(ORG[s], int(counts[s])) for s in species])

    rows = []
    for src in species:
        tr_mask = (df.organism == src).values
        tr_idx = np.where(tr_mask)[0]
        # hold out 20 % of the source species' clusters -> honest diagonal
        clusters = df.iloc[tr_idx].cluster.unique()
        rng = np.random.default_rng(SEED)
        holdout_c = set(rng.choice(clusters, size=max(2, int(0.2 * len(clusters))), replace=False))
        is_holdout = df.cluster.isin(holdout_c).values
        src_test_idx = tr_idx[is_holdout[tr_idx]]
        trainable = tr_idx[~is_holdout[tr_idx]]

        pos = df.iloc[trainable][df.iloc[trainable].label == 1]
        negs = df.iloc[trainable][df.iloc[trainable].label == 0]
        # training negatives: same-protein negatives of the source plus human background
        bg = df[(df.source == "background_cys") & (df.organism == "9606") & ~df.cluster.isin(holdout_c)]
        keep_bg = bg.sample(n=min(len(bg), len(pos)), random_state=SEED)
        fit_idx = np.concatenate([pos.index.values, negs.index.values, keep_bg.index.values])

        clusters_left = df.iloc[trainable].cluster.unique()
        rng = np.random.default_rng(SEED)
        val_c = set(rng.choice(clusters_left, size=max(2, int(0.15 * len(clusters_left))), replace=False))
        val_idx = trainable[df.iloc[trainable].cluster.isin(val_c).values]
        fit_idx = np.setdiff1d(fit_idx, val_idx)

        set_seed(SEED)
        rf = make_rf()
        rf.fit(np.hstack([Xflat[fit_idx], Cctx[fit_idx]]), df.label.values[fit_idx])
        net = EmbeddingMLP(in_dim=Xesm.shape[1])
        net, best, _ = train_torch(net, Xesm[fit_idx], Cctx[fit_idx],
                                   df.label.values[fit_idx].astype(np.float32),
                                   Xesm[val_idx], Cctx[val_idx],
                                   df.label.values[val_idx].astype(np.float32),
                                   epochs=30, lr=5e-4, patience=6)

        for tgt in species:
            te = src_test_idx if tgt == src else np.where((df.organism == tgt).values)[0]
            if len(te) < 40 or df.label.values[te].sum() < 20:
                continue
            p_rf = rf.predict_proba(np.hstack([Xflat[te], Cctx[te]]))[:, 1]
            net.eval()
            with torch.no_grad():
                p_esm = torch.sigmoid(net(torch.as_tensor(Xesm[te]),
                                          torch.as_tensor(Cctx[te]))).numpy()
            p = 1 / (1 + np.exp(-0.5 * (logit(p_rf) + logit(p_esm))))
            y = df.label.values[te]
            rows.append({"source": ORG[src], "target": ORG[tgt],
                         "auc": roc_auc_score(y, p), "n_pos": int(y.sum()),
                         "src_pos": int(len(pos)), "same_species": src == tgt})
        print(f"{ORG[src]:12s} done (train pos={len(pos)})", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "species_matrix.csv", index=False)
    mat = out.pivot(index="source", columns="target", values="auc")
    order = [ORG[s] for s in species if ORG[s] in mat.index]
    mat = mat.reindex(index=order, columns=order)
    print("\nAUC matrix (rows=train species, cols=test species)\n", mat.round(3).to_string())
    mat.to_csv(RESULTS / "species_matrix_auc.csv")
    (RESULTS / "species_matrix_summary.json").write_text(json.dumps(
        {"within_mean": float(np.nanmean(np.diag(mat.values))),
         "between_mean": float(np.nanmean(mat.values[~np.eye(len(mat), dtype=bool)])),
         "annotated": int((~np.isnan(mat.values)).sum())}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
