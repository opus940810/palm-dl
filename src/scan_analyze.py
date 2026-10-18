"""Analyse the proteome-wide scan: target-decoy FDR, shortlists, novelty.

Decoys are the same windows with the central cysteine mutated to alanine; the
random forest view is insensitive to the identity of the central residue
(§2.9), so decoys approximate the null distribution of a non-palmitoylatable
cysteine with identical flanking sequence.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from config import DATA_PROC, RESULTS

OUT = RESULTS / "scan"


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg style q-values on scores interpreted as p-values."""
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def main() -> None:
    chunks = sorted(OUT.glob("chunk_*.csv"))
    tgt = pd.concat([pd.read_csv(c) for c in chunks], ignore_index=True)
    meta = pd.read_csv(OUT / "targets.csv")
    tgt = tgt.merge(meta[["acc", "pos", "window", "compartment", "protein_length"]],
                    on=["acc", "pos"], how="left")
    dec = pd.read_csv(OUT / "decoys_shuffled.csv")
    print(f"targets {len(tgt)}, decoys {len(dec)}")

    # empirical target-decoy FDR at each target threshold
    tgt = tgt.sort_values("p", ascending=False).reset_index(drop=True)
    dec_sorted = np.sort(dec.p.values)[::-1]
    scale = len(tgt) / len(dec)
    n_dec_ge = np.searchsorted(-dec_sorted, -tgt.p.values, side="right")
    fdr = np.minimum(1.0, n_dec_ge * scale / np.arange(1, len(tgt) + 1))
    tgt["fdr"] = fdr
    tgt["q"] = np.minimum.accumulate(fdr[::-1])[::-1]

    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    known = set(zip(sites[(sites.label == 1) & sites.is_human].acc,
                    sites[(sites.label == 1) & sites.is_human].pos))
    lit = pd.read_csv(RESULTS / "literature" / "literature_sites.csv")
    lit_set = set(zip(lit.acc, lit.pos))
    tgt["known_uniprot"] = [(a, p) in known for a, p in zip(tgt.acc, tgt.pos)]
    tgt["in_literature"] = [(a, p) in lit_set for a, p in zip(tgt.acc, tgt.pos)]

    summary = {"n_targets": int(len(tgt)), "n_proteins": int(tgt.acc.nunique()),
               "n_known_uniprot": int(len(known)), "n_literature": int(len(lit_set))}
    rows = []
    for f in (0.01, 0.05, 0.10, 0.20, 0.30):
        sel = tgt[tgt.q <= f]
        known_hits = int(sel.known_uniprot.sum())
        lit_hits = int(sel.in_literature.sum())
        rows.append({"fdr": f, "n_predictions": int(len(sel)),
                     "n_proteins": int(sel.acc.nunique()),
                     "known_recovered": known_hits,
                     "known_recall": known_hits / max(len(known), 1),
                     "literature_recovered": lit_hits,
                     "novel": int(len(sel) - known_hits),
                     "threshold": float(sel.p.min()) if len(sel) else np.nan})
        if f in (0.01, 0.05, 0.20):
            sel.to_csv(OUT / f"shortlist_fdr{int(f * 100)}.csv", index=False)
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "fdr_summary.csv", index=False)
    print(tab.round(3).to_string(index=False))

    # enrichment of the 5% shortlist relative to all scanned cysteines
    sel = tgt[tgt.q <= 0.05]
    enrich = []
    if "compartment" in tgt.columns:
        for comp in tgt.compartment.unique():
            a = (sel.compartment == comp).mean() if "compartment" in sel else 0
            b = (tgt.compartment == comp).mean()
            enrich.append({"compartment": comp, "share_shortlist": a, "share_all": b,
                           "enrichment": a / b if b else np.nan})
        pd.DataFrame(enrich).to_csv(OUT / "compartment_enrichment.csv", index=False)
        print(pd.DataFrame(enrich).round(3).to_string(index=False))

    # per-protein multiplicity
    per_prot = sel.groupby("acc").size()
    summary.update({
        "shortlist_1pct": int((tgt.q <= 0.01).sum()),
        "shortlist_5pct": int(len(sel)),
        "proteins_5pct": int(sel.acc.nunique()),
        "median_sites_per_protein_5pct": float(per_prot.median()) if len(per_prot) else 0.0,
        "max_sites_per_protein_5pct": int(per_prot.max()) if len(per_prot) else 0,
    })
    (OUT / "scan_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    tgt.to_csv(OUT / "all_scores.csv", index=False)


if __name__ == "__main__":
    main()
