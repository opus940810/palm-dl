"""Curated subset of the literature-mined sites + overlap with the proteome scan.

Curation rules (applied by reading every sentence in review_sheet.tsv):
  accept  - the sentence explicitly attributes S-palmitoylation to that protein
            at that cysteine (including loss-of-palmitoylation mutants), and the
            mapped accessions match the protein named in the sentence;
  reject  - the site belongs to another protein in the sentence (enzyme token
            captured instead of substrate), a review table row with a mismatched
            protein, a catalytic-site mutant of a DHHC enzyme, a purely
            predicted site, a negative result, or a position never stated.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from config import RESULTS

OUT = RESULTS / "literature"
ACCEPT = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18, 19, 21, 22, 24, 25, 26, 28, 29, 30, 31,
          32, 33, 36, 37, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 53, 54, 55, 58, 59,
          60, 61, 63, 65, 69, 70, 71, 72, 75, 76, 77, 78, 79, 81, 82, 83, 84, 86, 87, 91, 92,
          93, 94, 95, 96, 97, 98, 101, 102, 103, 105, 107, 108, 109, 111, 112, 114, 115, 116,
          117, 119, 122, 123, 124, 125, 126, 127, 129, 130, 131, 132, 133, 134, 135]


def main() -> None:
    review = pd.read_csv(OUT / "review_sheet.tsv", sep="\t")
    review["curated"] = review.i.isin(ACCEPT)
    review[["i", "acc", "pos", "token", "protein", "pmcid", "curated"]].to_csv(
        OUT / "curated_sites.csv", index=False)
    cur = review[review.curated].copy()
    print(f"curated {len(cur)} / {len(review)} sites, {cur.acc.nunique()} proteins")

    # motif check on the curated subset (should beat the unfiltered 0.537)
    pred = pd.read_csv(OUT / "literature_predictions.csv")
    neg = pred[pred.label == 0]
    cur_pred = pred[(pred.label == 1)].merge(cur[["acc", "pos"]], on=["acc", "pos"], how="inner")
    if len(cur_pred) and len(neg):
        y = np.r_[np.ones(len(cur_pred)), np.zeros(len(neg))]
        auc_pwm = roc_auc_score(y, np.r_[cur_pred.pwm, neg.pwm])
        auc_model = roc_auc_score(y, np.r_[cur_pred.p, neg.p])
        print(f"curated vs same-protein negatives: PWM AUC {auc_pwm:.3f}, model AUC {auc_model:.3f}")

    # overlap with the proteome-wide shortlists
    out = {"n_curated": int(len(cur)), "n_proteins": int(cur.acc.nunique())}
    all_cys = pd.read_csv(RESULTS / "scan" / "all_scores.csv")
    total = len(all_cys)
    for tier in (1, 5, 20):
        path = RESULTS / "scan" / f"shortlist_fdr{tier}.csv"
        if not path.exists():
            continue
        short = pd.read_csv(path)
        keys = set(zip(short.acc, short.pos))
        hits = sum((a, p) in keys for a, p in zip(cur.acc, cur.pos))
        expected = len(cur) * len(short) / total
        out[f"fdr{tier}_hits"] = int(hits)
        out[f"fdr{tier}_expected"] = float(expected)
        out[f"fdr{tier}_enrichment"] = float(hits / expected) if expected else float("nan")
        print(f"FDR {tier}%: {hits} curated sites in shortlist of {len(short)} "
              f"(expected {expected:.2f}, enrichment {hits / max(expected, 1e-9):.1f}x)")
    (OUT / "curated_overlap.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
