"""Shuffled-flank decoys for proteome-wide FDR, scored with the forest view.

Mutating the centre cysteine (C->A) does not destroy the flanking motif, and the
model is deliberately insensitive to the identity of the centre residue because
every candidate is a cysteine. A null that asks "what would a cysteine with no
motif score?" must therefore shuffle the 30 flanking residues while keeping the
centre cysteine in place.

The scan uses the forest view (the ESM view could not be run at proteome scale on
this CPU host); probabilities come from the forest Platt calibrator fitted on the
same out-of-fold validation rows.
"""
from __future__ import annotations

import joblib
import numpy as np
import pandas as pd

from config import MODELS, RESULTS
from features import context_features, encode_flat
from final_model import logit, load

OUT = RESULTS / "scan"


def main() -> None:
    chunks = sorted(OUT.glob("chunk_*.csv"))
    tgt = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True)
    meta = pd.read_csv(OUT / "targets.csv")
    tgt = tgt.merge(meta[["acc", "pos", "window", "compartment", "protein_length"]],
                    on=["acc", "pos"], how="left")
    rng = np.random.default_rng(1)
    idx = rng.choice(len(tgt), size=20000, replace=False)
    d = tgt.iloc[idx].copy()

    def shuffle_flanks(w: str) -> str:
        chars = list(w[:15] + w[16:])
        rng.shuffle(chars)
        return "".join(chars[:15]) + w[15] + "".join(chars[15:])

    d["window"] = [shuffle_flanks(w) for w in d.window]
    d["compartment"] = d["compartment"].fillna("unknown")
    d["protein_length"] = d["protein_length"].fillna(500)

    rf, net, cal, meta_model = load()
    cal_rf = joblib.load(MODELS / "calibrator_rf.joblib")
    X = encode_flat(d.window.tolist())
    C = context_features(d)
    p_rf = rf.predict_proba(np.hstack([X, C]))[:, 1]
    p = cal_rf.predict_proba(logit(p_rf).reshape(-1, 1))[:, 1]
    pd.DataFrame({"acc": d.acc.values, "pos": d.pos.values, "p": p}).to_csv(
        OUT / "decoys_shuffled.csv", index=False)
    print("shuffled decoys:", len(d), "median p", round(float(np.median(p)), 4),
          "p99", round(float(np.quantile(p, 0.99)), 3))


if __name__ == "__main__":
    main()
