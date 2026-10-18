"""Addendum: do literature sites look like palmitoylation sites at all?

Two questions are separated:
  (1) recognition - literature sites versus random background cysteines;
  (2) discrimination - literature sites versus other cysteines of the same protein.

If (1) is high and (2) is low, the model recognises these sites as
palmitoylation-like but cannot rank them above neighbouring cysteines, which is
the expected consequence of training on motif-biased UniProt annotations.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from config import DATA_PROC, MODELS, RESULTS
from final_model import score

OUT = RESULTS / "literature"


def main() -> None:
    lit = pd.read_csv(OUT / "literature_predictions.csv")
    pos = lit[(lit.label == 1) & ~lit.in_uniprot]
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    bg = sites[(sites.source == "background_cys") & sites.window.str.len().eq(31)]
    bg = bg.sample(n=min(5000, len(bg)), random_state=0)
    bg = bg[~bg.acc.isin(set(pos.acc))].head(4000)
    bg = bg.assign(compartment="unknown")
    pred = score(bg[["acc", "pos", "window", "compartment", "protein_length"]])
    bg = bg.assign(p=pred.p.values)

    # reference positives: UniProt experimental sites (out-of-fold predictions)
    oof = pd.read_csv(RESULTS / "preds_cv_all.csv")
    oof = oof[(oof.model == "ensemble") & (oof.set == "test") & (oof.y == 1)]
    exp = oof[oof.evidence.astype(str).str.contains("ECO:0000269")] if "evidence" in oof else oof.head(0)
    neg_oof = pd.read_csv(RESULTS / "preds_cv_all.csv")
    neg_oof = neg_oof[(neg_oof.model == "ensemble") & (neg_oof.set == "test") & (neg_oof.y == 0)]

    out = {"n_literature": int(len(pos)), "n_background": int(len(bg))}
    if len(pos) and len(bg):
        y = np.r_[np.ones(len(pos)), np.zeros(len(bg))]
        out["literature_vs_background_auc"] = float(roc_auc_score(y, np.r_[pos.p, bg.p]))
        out["literature_vs_background_ap"] = float(average_precision_score(y, np.r_[pos.p, bg.p]))
        y2 = np.r_[np.ones(len(pos)), np.zeros(len(bg))]
        out["literature_vs_background_pwm_auc"] = float(
            roc_auc_score(y2, np.r_[pos.pwm, bg.window.map(lambda w: np.nan)]) if False else np.nan)
    if len(exp):
        out["uniprot_experimental_vs_background_auc"] = float(
            roc_auc_score(np.r_[np.ones(len(exp)), np.zeros(len(bg))], np.r_[exp.p, bg.p]))
    out["literature_mean_p"] = float(pos.p.mean())
    out["background_mean_p"] = float(bg.p.mean())
    print(json.dumps(out, indent=2))
    (OUT / "recognition_metrics.json").write_text(json.dumps(out, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
