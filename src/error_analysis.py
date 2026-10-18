"""Error structure and FDR-controlled operating points.

1. Discrimination by evidence tier and by annotated topology.
2. Operating points at fixed estimated false-discovery rate (FDR), derived from
   the inner validation folds and then measured on the outer test folds.
3. In-silico mutagenesis of the random-forest view: which positions and which
   substitutions change the predicted probability of a true site.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from config import DATA_PROC, RESULTS, WINDOW
from features import context_features, encode_flat
from final_model import load

OUT = RESULTS / "error_analysis"
OUT.mkdir(exist_ok=True)


def tier_auc(df: pd.DataFrame, exp_mask: pd.Series) -> dict:
    pos = df[df.y == 1]
    neg = df[df.y == 0]
    out = {}
    for name, flag in (("experimental", True), ("inferred", False)):
        sub = pd.concat([pos[exp_mask.loc[pos.index].values == flag], neg])
        if sub.y.nunique() < 2:
            continue
        out[name] = {"n_pos": int(sub.y.sum()), "auc": roc_auc_score(sub.y, sub.p),
                     "aupr": average_precision_score(sub.y, sub.p)}
    return out


def fdr_operating_point(va: pd.DataFrame, te: pd.DataFrame, target_fdr: float) -> dict:
    """Threshold on validation whose estimated FDR is <= target, then test metrics."""
    y, p = va.y.values, va.p.values
    order = np.argsort(-p)
    y = y[order]
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    fdr = fp / np.maximum(tp + fp, 1)
    ok = np.where(fdr <= target_fdr)[0]
    if len(ok) == 0:
        return {}
    k = ok.max()
    thr = np.sort(p)[::-1][k]
    pred = te.p.values >= thr
    tp_t = int(((te.y == 1) & pred).sum())
    fp_t = int(((te.y == 0) & pred).sum())
    return {"target_fdr": target_fdr, "threshold": float(thr),
            "val_fdr": float(fdr[k]), "test_precision": tp_t / max(tp_t + fp_t, 1),
            "test_recall": tp_t / int((te.y == 1).sum()), "n_predicted": int(pred.sum())}


def mutagenesis(n_sites: int = 300, seed: int = 0) -> None:
    rf, _, _, _ = load()
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    pos = sites[(sites.label == 1) & sites.window.str.len().eq(WINDOW) & sites.is_human]
    pos = pos.sample(n=min(n_sites, len(pos)), random_state=seed)
    windows = pos.window.tolist()
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    ctx_row = context_features(pd.DataFrame({"compartment": ["unknown"] * len(windows),
                                             "protein_length": [500] * len(windows),
                                             "pos": [250] * len(windows)}))[:1]
    base_X = encode_flat(windows)
    base = rf.predict_proba(np.hstack([base_X, np.repeat(ctx_row, len(windows), axis=0)]))[:, 1]
    eff = np.zeros((20, WINDOW))
    for j in range(WINDOW):
        mut = np.array([list(w) for w in windows])
        for a_i, a in enumerate(aa):
            mut[:, j] = a
            X = encode_flat(["".join(row) for row in mut])
            p = rf.predict_proba(np.hstack([X, np.repeat(ctx_row, len(windows), axis=0)]))[:, 1]
            eff[a_i, j] = float(np.mean(p - base))
    np.save(OUT / "mutagenesis_effect.npy", eff)
    pd.DataFrame(eff, index=aa, columns=[str(j - WINDOW // 2) for j in range(WINDOW)]).to_csv(
        OUT / "mutagenesis_effect.csv")
    print("mutagenesis done; strongest effects:")
    flat = [(eff[i, j], aa[i], j - WINDOW // 2) for i in range(20) for j in range(WINDOW)]
    for v, a, j in sorted(flat, reverse=True)[:8]:
        print(f"  {a} at {j:+d}: {v:+.3f}")
    for v, a, j in sorted(flat)[:8]:
        print(f"  {a} at {j:+d}: {v:+.3f}")


def main() -> None:
    d = pd.read_csv(RESULTS / "preds_cv_all.csv")
    d = d[(d.model == "ensemble") & (d.set == "test")]
    va = pd.read_csv(RESULTS / "preds_cv_all.csv")
    va = va[(va.model == "ensemble") & (va.set == "val")]
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    meta_cols = sites[["acc", "pos", "evidence", "compartment", "source"]]
    d = d.merge(meta_cols, on=["acc", "pos"], how="left", suffixes=("", "_s"))
    va = va.merge(meta_cols, on=["acc", "pos"], how="left", suffixes=("", "_s"))

    exp = d.evidence.astype(str).str.contains("ECO:0000269")
    res = {"tier": tier_auc(d, exp)}
    print("tier:", json.dumps(res["tier"], indent=2))

    comp_rows = []
    for comp in ["cytoplasmic", "extracellular", "lumenal", "unknown"]:
        sub = d[d.compartment == comp]
        if sub.y.nunique() > 1 and sub.y.sum() >= 10:
            comp_rows.append({"compartment": comp, "n": len(sub), "n_pos": int(sub.y.sum()),
                              "auc": roc_auc_score(sub.y, sub.p),
                              "aupr": average_precision_score(sub.y, sub.p),
                              "mean_p_pos": float(sub[sub.y == 1].p.mean()),
                              "mean_p_neg": float(sub[sub.y == 0].p.mean())})
    pd.DataFrame(comp_rows).to_csv(OUT / "compartment_metrics.csv", index=False)
    print(pd.DataFrame(comp_rows).round(3).to_string(index=False))

    fdr_rows = [fdr_operating_point(va, d, f) for f in (0.01, 0.05, 0.10, 0.20)]
    fdr = pd.DataFrame([r for r in fdr_rows if r]).sort_values("target_fdr")
    fdr.to_csv(OUT / "fdr_operating_points.csv", index=False)
    print("\nFDR operating points:\n", fdr.round(3).to_string(index=False))

    mutagenesis()
    (OUT / "summary.json").write_text(json.dumps({"tier": res["tier"],
                                                  "fdr": fdr.to_dict(orient="records")}, indent=2),
                                      encoding="utf-8")


if __name__ == "__main__":
    main()
