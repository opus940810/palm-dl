"""Cross-species transfer matrix v2: frozen target sets, matched negatives.

For each source species one model is trained on that species' positives plus its
same-protein negatives; human background negatives are excluded so that every
row uses the same negative class. Twenty percent of the source clusters are held
out for the diagonal, and another 15% for early stopping. Target rows are
restricted to positives plus same-protein negatives and are removed if their
cluster appears anywhere in training (fit or validation), so no target family is
seen during fitting or early stopping. Off-diagonal cells use the same model as
the diagonal (a single model per source), which the manuscript states explicitly.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, EMB, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from final_model import logit
from models import EmbeddingMLP, make_rf, set_seed, train_torch

torch.set_num_threads(int(os.environ.get("PALMDL_THREADS", "1")))
ORG = {"9606": "Human", "10090": "Mouse", "10116": "Rat", "9913": "Cattle",
       "3702": "A. thaliana", "9823": "Pig"}
MIN_POS = 85
MIN_TEST_POS = 20


def main() -> None:
    set_seed(SEED)
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")
    smap = sites.set_index(["acc", "pos"]).source
    df = split.merge(sites[["acc", "pos", "window", "compartment", "protein_length"]],
                     on=["acc", "pos"], validate="one_to_one")
    df["source"] = [smap.get((a, p)) for a, p in zip(df.acc, df.pos)]
    df["organism"] = df["organism"].astype(str).str.split(".").str[0]
    df = df[df.window.str.len() == WINDOW].reset_index(drop=True)
    Xflat = encode_flat(df.window.tolist())
    Cctx = context_features(df)
    Xesm = np.hstack([np.load(EMB / "esm2_mean.npy").astype(np.float32),
                      np.load(EMB / "esm2_center.npy").astype(np.float32)])

    counts = df[df.label == 1].organism.value_counts()
    species = [s for s in ORG if s in counts.index and counts[s] >= MIN_POS]
    print("species:", [(ORG[s], int(counts[s])) for s in species], flush=True)

    rows = []
    for src in species:
        src_rows = df[df.organism == src]
        src_clusters = src_rows.cluster.unique()
        rng = np.random.default_rng(SEED)
        holdout = set(rng.choice(src_clusters, size=max(2, int(0.2 * len(src_clusters))), replace=False))
        left = [c for c in src_clusters if c not in holdout]
        val_c = set(rng.choice(left, size=max(2, int(0.15 * len(left))), replace=False))
        fit_mask = (df.organism == src) & (~df.cluster.isin(holdout)) & (~df.cluster.isin(val_c))
        fit_mask &= df.source != "background_cys"
        val_mask = (df.organism == src) & df.cluster.isin(val_c)
        fit_idx = np.where(fit_mask.values)[0]
        val_idx = np.where(val_mask.values)[0]
        used_clusters = set(df.iloc[fit_idx].cluster.unique()) | set(df.iloc[val_idx].cluster.unique())
        print(f"{ORG[src]}: fit={len(fit_idx)} (pos={int(df.iloc[fit_idx].label.sum())}) "
              f"val={len(val_idx)}", flush=True)

        set_seed(SEED)
        rf = make_rf()
        rf.fit(np.hstack([Xflat[fit_idx], Cctx[fit_idx]]), df.label.values[fit_idx])
        net = EmbeddingMLP(in_dim=Xesm.shape[1])
        net, best, _ = train_torch(net, Xesm[fit_idx], Cctx[fit_idx],
                                   df.label.values[fit_idx].astype(np.float32),
                                   Xesm[val_idx], Cctx[val_idx],
                                   df.label.values[val_idx].astype(np.float32),
                                   epochs=30, lr=5e-4, patience=6)
        print(f"  esm stopped epoch {best['epoch']} val_auc {best['score']:.4f}", flush=True)

        for tgt in species:
            cand = df[(df.organism == tgt) & (df.source != "background_cys")]
            keep = ~cand.cluster.isin(used_clusters)
            te_idx = cand.index[keep].values
            if tgt == src:
                te_idx = df.index[(df.organism == src) & df.cluster.isin(holdout)
                                  & (df.source != "background_cys")].values
            removed = int(keep.sum() == 0 and tgt != src)
            y = df.label.values[te_idx]
            if len(te_idx) < 30 or y.sum() < MIN_TEST_POS:
                print(f"  {ORG[src]}->{ORG[tgt]}: skipped (n={len(te_idx)}, pos={int(y.sum())})", flush=True)
                rows.append({"source": ORG[src], "target": ORG[tgt], "auc": np.nan,
                             "n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum()),
                             "n_removed": int((~keep).sum()) if tgt != src else 0,
                             "note": "insufficient test positives"})
                continue
            p_rf = rf.predict_proba(np.hstack([Xflat[te_idx], Cctx[te_idx]]))[:, 1]
            net.eval()
            with torch.no_grad():
                p_esm = torch.sigmoid(net(torch.as_tensor(Xesm[te_idx]),
                                          torch.as_tensor(Cctx[te_idx]))).numpy()
            p = 1 / (1 + np.exp(-0.5 * (logit(p_rf) + logit(p_esm))))
            auc = roc_auc_score(y, p)
            rows.append({"source": ORG[src], "target": ORG[tgt], "auc": auc,
                         "n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum()),
                         "n_removed": int((~keep).sum()) if tgt != src else 0, "note": ""})
            print(f"  {ORG[src]}->{ORG[tgt]}: AUC {auc:.3f} "
                  f"(pos={int(y.sum())}, neg={int(len(y) - y.sum())})", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "species_matrix_v2.csv", index=False)
    mat = out.pivot(index="source", columns="target", values="auc")
    order = [ORG[s] for s in species if ORG[s] in mat.index]
    mat = mat.reindex(index=order, columns=order)
    mat.to_csv(RESULTS / "species_matrix_v2_auc.csv")
    print("\n", mat.round(3).to_string(), flush=True)
    off = mat.values[~np.eye(len(mat), dtype=bool)]
    (RESULTS / "species_matrix_v2_summary.json").write_text(json.dumps({
        "within_mean": float(np.nanmean(np.diag(mat.values))),
        "between_mean": float(np.nanmean(off)),
        "between_min": float(np.nanmin(off)),
        "between_max": float(np.nanmax(off)),
        "annotated": int((~np.isnan(mat.values)).sum()),
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
