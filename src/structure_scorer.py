"""Structure-aware PalmDL prototype.

Design (chosen after three comparisons):
  * RF view with structure:     0.890 -> 0.910  (+0.020, 5/5 folds)
  * reranking the full ensemble: 0.919 -> 0.914  (fails: the ensemble already
    carries most of the structural signal, and a 4-parameter reranker overfits)
  * structure as a third view:  0.918 -> 0.924  (+0.006, 5/5 folds)  <- adopted

So the prototype replaces the RF view with an RF trained on one-hot + context +
AlphaFold features, keeps the frozen ESM view, and averages the two logits.
It can only score proteins for which an AlphaFold model exists.

Stages
  fit         train the structure-aware RF on all structured human sites
  literature  apply to the curated literature sites (domain-shift test)
  scan        apply to proteome-scan cysteines with an AlphaFold model
"""
from __future__ import annotations

import argparse
import json

import joblib
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, MODELS, RESULTS, SEED, WINDOW
from features import context_features, encode_flat
from final_model import embed_windows, load, logit
from models import make_rf, set_seed

OUT = RESULTS / "structure"
FEATS = ["plddt_site", "plddt_window", "packing"]
RF_STRUCT = MODELS / "rf_struct.joblib"


def fit() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    st = pd.read_csv(OUT / "site_structure_features.csv")
    d = sites.merge(st, on=["acc", "pos"], how="inner")
    X = np.hstack([encode_flat(d.window.tolist()), context_features(d), d[FEATS].values])
    set_seed(SEED)
    rf = make_rf().fit(X, d.label.values)
    joblib.dump(rf, RF_STRUCT)
    meta = {"n_fit": int(len(d)), "n_pos": int(d.label.sum()), "features": FEATS,
            "n_proteins": int(d.acc.nunique())}
    (MODELS / "rf_struct_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


def score_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Score rows that already carry structure features."""
    rf_base, net, cal, _ = load()
    rf_struct = joblib.load(RF_STRUCT)
    X = encode_flat(df.window.tolist())
    C = context_features(df)
    S = df[FEATS].values.astype(float)
    p_rf = rf_base.predict_proba(np.hstack([X, C]))[:, 1]
    p_rs = rf_struct.predict_proba(np.hstack([X, C, S]))[:, 1]
    E = embed_windows(df.window.tolist(), bs=192)
    with torch.no_grad():
        p_esm = torch.sigmoid(net(torch.as_tensor(E), torch.as_tensor(C))).numpy()
    out = df.copy()
    out["p_rf"] = p_rf
    out["p_rf_struct"] = p_rs
    out["p_esm"] = p_esm
    out["p_seq"] = 1 / (1 + np.exp(-0.5 * (logit(p_rf) + logit(p_esm))))
    out["p_struct"] = 1 / (1 + np.exp(-0.5 * (logit(p_rs) + logit(p_esm))))
    return out


def literature() -> None:
    from concurrent.futures import ThreadPoolExecutor

    from af_features import fetch, features as af_features
    from build_dataset import compartment as comp_of
    from build_dataset import parse_topology
    from build_dataset import window as mkwin

    lit = pd.read_csv(RESULTS / "literature" / "literature_predictions.csv")
    lit = lit[["acc", "pos", "label"]].copy()
    seqs, topo = {}, {}
    for name in ("uniprot_negatives.tsv", "uniprot_lipid_reviewed.tsv"):
        raw = pd.read_csv(DATA_PROC.parent / "raw" / name, sep="\t", dtype=str).fillna("")
        for acc, s, tp in zip(raw["Entry"], raw["Sequence"], raw.get("Topological domain", [""] * len(raw))):
            seqs[acc] = s.upper()
            topo[acc] = tp
    rows = []
    for _, r in lit.iterrows():
        s = seqs.get(r.acc)
        if not s or r.pos < 1 or r.pos > len(s):
            continue
        rows.append({"acc": r.acc, "pos": int(r.pos), "label": int(r.label),
                     "window": mkwin(s, int(r.pos)), "protein_length": len(s),
                     "compartment": comp_of(int(r.pos), parse_topology(topo.get(r.acc, "")), None)})
    d0 = pd.DataFrame(rows)
    with ThreadPoolExecutor(max_workers=6) as ex:
        list(ex.map(fetch, sorted(d0.acc.unique())))
    have = {p.stem for p in (OUT / "af").glob("*.pdb")}
    d0 = d0[d0.acc.isin(have)]
    af = af_features(d0)
    if af.empty:
        print("no AlphaFold models for literature sites")
        return
    d = af.merge(d0[["acc", "pos", "window", "protein_length", "compartment", "label"]],
                 on=["acc", "pos"], how="left")
    print(f"literature rows with structure: {len(d)} ({int(d.label.sum())} positives)")
    s = score_frame(d)
    s.to_csv(OUT / "literature_struct_scored.csv", index=False)
    res = {}
    if s.label.nunique() == 2:
        res = {"auc_seq": float(roc_auc_score(s.label, s.p_seq)),
               "auc_struct": float(roc_auc_score(s.label, s.p_struct)),
               "delta": float(roc_auc_score(s.label, s.p_struct) - roc_auc_score(s.label, s.p_seq)),
               "n_pos": int(s.label.sum()), "n": int(len(s))}
        print(json.dumps(res, indent=2))
    (OUT / "literature_struct_eval.json").write_text(json.dumps(res, indent=2), encoding="utf-8")


def scan() -> None:
    from af_features import features as af_features

    have = {p.stem for p in (OUT / "af").glob("*.pdb")}
    tgt = pd.read_csv(RESULTS / "scan" / "targets.csv")
    tgt = tgt[tgt.acc.isin(have)]
    print(f"scan targets with structure: {len(tgt)}")
    feats = af_features(tgt[["acc", "pos", "window", "protein_length", "compartment"]])
    meta = tgt[["acc", "pos", "window", "protein_length", "compartment"]]
    scores = pd.read_csv(RESULTS / "scan" / "all_scores.csv")[["acc", "pos", "known_uniprot"]]
    d = feats.merge(meta, on=["acc", "pos"], how="left").merge(scores, on=["acc", "pos"], how="left")
    s = score_frame(d)
    s.to_csv(OUT / "scan_struct_scored.csv", index=False)
    res = {"n_sites": int(len(s)), "n_known": int(s.known_uniprot.fillna(False).sum())}
    for k in (322, 668):
        a = set(s.nlargest(k, "p_seq")[["acc", "pos"]].itertuples(index=False, name=None))
        b = set(s.nlargest(k, "p_struct")[["acc", "pos"]].itertuples(index=False, name=None))
        known = dict(zip(zip(s.acc, s.pos), s.known_uniprot.fillna(False)))
        res[f"top{k}_known_seq"] = int(sum(known[x] for x in a))
        res[f"top{k}_known_struct"] = int(sum(known[x] for x in b))
        res[f"top{k}_overlap"] = int(len(a & b))
        print(f"top {k}: known seq {res[f'top{k}_known_seq']} -> struct {res[f'top{k}_known_struct']}; "
              f"overlap {res[f'top{k}_overlap']}")
    (OUT / "scan_struct_eval.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["fit", "literature", "scan", "all"])
    args = ap.parse_args()
    if args.stage in ("fit", "all"):
        fit()
    if args.stage in ("literature", "all"):
        literature()
    if args.stage in ("scan", "all"):
        scan()


if __name__ == "__main__":
    main()
