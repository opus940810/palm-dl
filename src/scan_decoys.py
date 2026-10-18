"""Better decoy set for proteome-wide FDR: shuffled flanks, centre cysteine kept.

Mutating the centre cysteine (C->A) does not destroy the flanking motif, and the
model is deliberately insensitive to the identity of the centre residue because
every candidate is a cysteine.  A null that asks "what would a cysteine with no
motif score?" must therefore shuffle the 30 flanking residues while keeping the
centre cysteine in place.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import RESULTS
from final_model import score

OUT = RESULTS / "scan"


def main() -> None:
    tgt = pd.read_csv(OUT / "all_scores.csv")
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
    pred = score(d[["acc", "pos", "window", "compartment", "protein_length"]])
    pd.DataFrame({"acc": d.acc.values, "pos": d.pos.values, "p": pred.p.values}).to_csv(
        OUT / "decoys_shuffled.csv", index=False)
    print("shuffled decoys:", len(d), "median p", round(float(pred.p.median()), 4),
          "p99", round(float(pred.p.quantile(0.99)), 3))


if __name__ == "__main__":
    main()
