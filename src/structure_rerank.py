"""Structure-aware reranker prototype.

The deployed PalmDL score is sequence-only.  Here we learn a small logistic
reranker on top of it using three AlphaFold features (pLDDT at the cysteine,
mean pLDDT over +/-5 residues, C-alpha packing density), trained and evaluated
with the existing cluster-disjoint folds.

Stages
  cv        paired CV of the reranker against the sequence score
  fit       fit and save the reranker on all structured human sites
  literature  apply to the curated literature sites (domain-shift test)
  scan      apply to the proteome scan for proteins that have an AlphaFold model
"""
from __future__ import annotations

import argparse
import json

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from config import DATA_PROC, MODELS, RESULTS, SEED

OUT = RESULTS / "structure"
FEATS = ["plddt_site", "plddt_window", "packing"]


def logit(p, eps=1e-6):
    p = np.clip(np.asarray(p, dtype=float), eps, 1 - eps)
    return np.log(p / (1 - p))


def structured_sites() -> pd.DataFrame:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    split = pd.read_csv(DATA_PROC / "splits_cv.csv")[["acc", "pos", "fold"]]
    st = pd.read_csv(OUT / "site_structure_features.csv")
    d = split.merge(sites[["acc", "pos", "label"]], on=["acc", "pos"])
    d = d.merge(st, on=["acc", "pos"], how="inner")
    oof = pd.read_csv(RESULTS / "preds_cv_all.csv")
    oof = oof[(oof.model == "ensemble") & (oof.set == "test")][["acc", "pos", "p", "y"]]
    d = d.merge(oof, on=["acc", "pos"], how="inner", suffixes=("", "_oof"))
    return d


def design(d: pd.DataFrame) -> np.ndarray:
    X = [logit(d.p.values)] + [d[f].values.astype(float) for f in FEATS]
    return np.column_stack(X)


def cv() -> dict:
    d = structured_sites()
    print(f"structured sites with OOF scores: {len(d)} ({int(d.label.sum())} positives, "
          f"{d.acc.nunique()} proteins)")
    rows = []
    for k in sorted(d.fold.unique()):
        tr = d[d.fold != k]
        te = d[d.fold == k]
        if te.label.nunique() < 2:
            continue
        lr = LogisticRegression(C=1.0, max_iter=2000).fit(design(tr), tr.label.values)
        p_re = lr.predict_proba(design(te))[:, 1]
        rows.append({"fold": k, "n": len(te), "n_pos": int(te.label.sum()),
                     "auc_seq": roc_auc_score(te.label.values, te.p.values),
                     "auc_rerank": roc_auc_score(te.label.values, p_re)})
        print(f"fold {k}: seq {rows[-1]['auc_seq']:.3f} -> rerank {rows[-1]['auc_rerank']:.3f}")
    res = pd.DataFrame(rows)
    delta = (res.auc_rerank - res.auc_seq)
    summary = {"n_sites": int(len(d)), "n_pos": int(d.label.sum()),
               "auc_seq": float(res.auc_seq.mean()), "auc_rerank": float(res.auc_rerank.mean()),
               "mean_delta": float(delta.mean()), "positive_folds": int((delta > 0).sum()),
               "folds": int(len(res))}
    res.to_csv(OUT / "reranker_cv.csv", index=False)
    (OUT / "reranker_cv_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


def fit() -> None:
    d = structured_sites()
    lr = LogisticRegression(C=1.0, max_iter=2000).fit(design(d), d.label.values)
    joblib.dump({"model": lr, "features": ["logit_p_seq"] + FEATS}, MODELS / "struct_reranker.joblib")
    meta = {"n_fit": int(len(d)), "n_pos": int(d.label.sum()),
            "plddt_pos_mean": float(d[d.label == 1].plddt_site.mean()),
            "plddt_neg_mean": float(d[d.label == 0].plddt_site.mean()),
            "packing_pos_mean": float(d[d.label == 1].packing.mean()),
            "packing_neg_mean": float(d[d.label == 0].packing.mean())}
    (MODELS / "struct_reranker_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


def literature() -> None:
    from af_features import features as af_features

    lit = pd.read_csv(RESULTS / "literature" / "literature_predictions.csv")
    lit = lit[["acc", "pos", "label", "p", "pwm"]].copy()
    af = af_features(lit[["acc", "pos"]].assign(
        window="X", protein_length=500, compartment="unknown"))
    if af.empty:
        print("no AlphaFold models for literature sites")
        return
    d = lit.merge(af, on=["acc", "pos"], how="inner")
    print(f"literature sites with structure: {len(d)} ({int(d.label.sum())} positives)")
    rr = joblib.load(MODELS / "struct_reranker.joblib")["model"]
    X = np.column_stack([logit(d.p.values)] + [d[f].values.astype(float) for f in FEATS])
    d["p_rerank"] = rr.predict_proba(X)[:, 1]
    out = {}
    if d.label.nunique() == 2:
        out["auc_seq"] = float(roc_auc_score(d.label, d.p))
        out["auc_rerank"] = float(roc_auc_score(d.label, d.p_rerank))
        print(f"literature: seq AUC {out['auc_seq']:.3f} -> rerank {out['auc_rerank']:.3f}")
    (OUT / "reranker_literature.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    d.to_csv(OUT / "reranker_literature_predictions.csv", index=False)


def scan() -> None:
    from pathlib import Path

    af_dir = OUT / "af"
    have = {p.stem for p in af_dir.glob("*.pdb")}
    tgt = pd.read_csv(RESULTS / "scan" / "targets.csv")
    tgt = tgt[tgt.acc.isin(have)].copy()
    print(f"scan targets with AlphaFold models: {len(tgt)} in {tgt.acc.nunique()} proteins")
    if tgt.empty:
        return
    from af_features import features as af_features

    feats = af_features(tgt[["acc", "pos", "window", "protein_length", "compartment"]])
    scores = pd.read_csv(RESULTS / "scan" / "all_scores.csv")
    d = feats.merge(scores[["acc", "pos", "p", "known_uniprot"]], on=["acc", "pos"], how="inner")
    rr = joblib.load(MODELS / "struct_reranker.joblib")["model"]
    X = np.column_stack([logit(d.p.values)] + [d[f].values.astype(float) for f in FEATS])
    d["p_rerank"] = rr.predict_proba(X)[:, 1]
    d["delta"] = d.p_rerank - d.p

    out = {"n_sites": int(len(d)), "n_known": int(d.known_uniprot.sum())}
    for k in (322, 668):
        top_seq = set(d.nlargest(k, "p")[["acc", "pos"]].itertuples(index=False, name=None))
        top_re = set(d.nlargest(k, "p_rerank")[["acc", "pos"]].itertuples(index=False, name=None))
        known_seq = sum(d[(d.acc == a) & (d.pos == p)].known_uniprot.sum() for a, p in top_seq)
        known_re = sum(d[(d.acc == a) & (d.pos == p)].known_uniprot.sum() for a, p in top_re)
        out[f"top{k}_known_seq"] = int(known_seq)
        out[f"top{k}_known_rerank"] = int(known_re)
        out[f"top{k}_overlap"] = int(len(top_seq & top_re))
        print(f"top {k}: known {known_seq} -> {known_re}; overlap {len(top_seq & top_re)}")
    out["mean_delta_novel"] = float(d[~d.known_uniprot].delta.mean())
    d.to_csv(OUT / "reranker_scan_predictions.csv", index=False)
    (OUT / "reranker_scan.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


def thirdview() -> None:
    """Is structure useful as a third view on top of the RF+ESM ensemble?"""
    from features import context_features, encode_flat
    from models import make_rf, set_seed

    d = structured_sites()
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    meta = sites[["acc", "pos", "window", "compartment", "protein_length"]]
    d = d.merge(meta, on=["acc", "pos"], how="left")
    esm = pd.read_csv(RESULTS / "preds_cv_esm.csv")
    esm = esm[(esm.model == "esm") & (esm.set == "test")][["acc", "pos", "p"]]
    d = d.merge(esm.rename(columns={"p": "p_esm"}), on=["acc", "pos"], how="inner")
    print(f"third-view rows: {len(d)} ({int(d.label.sum())} positives)")

    X = encode_flat(d.window.tolist())
    C = context_features(d)
    S = d[FEATS].values.astype(float)
    S = (S - S.mean(0)) / (S.std(0) + 1e-6)
    rows = []
    for k in sorted(d.fold.unique()):
        tr = np.where(d.fold.values != k)[0]
        te = np.where(d.fold.values == k)[0]
        if d.label.values[te].min() == d.label.values[te].max():
            continue
        set_seed(SEED + k)
        rf = make_rf().fit(np.hstack([X[tr], C[tr]]), d.label.values[tr])
        p_rf = rf.predict_proba(np.hstack([X[te], C[te]]))[:, 1]
        set_seed(SEED + k)
        rfs = make_rf().fit(np.hstack([X[tr], C[tr], S[tr]]), d.label.values[tr])
        p_rfs = rfs.predict_proba(np.hstack([X[te], C[te], S[te]]))[:, 1]
        p_esm = d.p_esm.values[te]
        y = d.label.values[te]
        seq = 1 / (1 + np.exp(-0.5 * (logit(p_rf) + logit(p_esm))))
        struct = 1 / (1 + np.exp(-0.5 * (logit(p_rfs) + logit(p_esm))))
        rows.append({"fold": k, "n": len(te), "auc_seq_ensemble": roc_auc_score(y, seq),
                     "auc_struct_ensemble": roc_auc_score(y, struct),
                     "auc_rf": roc_auc_score(y, p_rf), "auc_rf_struct": roc_auc_score(y, p_rfs)})
        print(f"fold {k}: RF {rows[-1]['auc_rf']:.3f}->{rows[-1]['auc_rf_struct']:.3f} | "
              f"ensemble {rows[-1]['auc_seq_ensemble']:.3f}->{rows[-1]['auc_struct_ensemble']:.3f}")
    res = pd.DataFrame(rows)
    delta = res.auc_struct_ensemble - res.auc_seq_ensemble
    summary = {"auc_seq_ensemble": float(res.auc_seq_ensemble.mean()),
               "auc_struct_ensemble": float(res.auc_struct_ensemble.mean()),
               "delta": float(delta.mean()), "positive_folds": int((delta > 0).sum()),
               "folds": int(len(res))}
    res.to_csv(OUT / "thirdview_cv.csv", index=False)
    (OUT / "thirdview_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["cv", "fit", "thirdview", "literature", "scan", "all"])
    args = ap.parse_args()
    if args.stage in ("cv", "all"):
        cv()
    if args.stage in ("fit", "all"):
        fit()
    if args.stage in ("thirdview", "all"):
        thirdview()
    if args.stage in ("literature", "all"):
        literature()
    if args.stage in ("scan", "all"):
        scan()


if __name__ == "__main__":
    main()
