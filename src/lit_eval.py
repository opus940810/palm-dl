"""Stage 3-4: resolve literature mentions to UniProt, verify cysteines, score.

Protein resolution uses the locally cached UniProt tables (23,897 reviewed
entries with gene names and sequences) instead of the REST API, which is
intermittent.  Positives are literature-mined sites whose residue at the
mentioned position is a cysteine in the reviewed sequence.  Negatives are other
cysteines of the same protein that are neither mentioned nor annotated as
S-palmitoylated in UniProt 2026_03.
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from build_dataset import window
from config import AA, DATA_PROC, DATA_RAW, MODELS, RESULTS, WINDOW
from final_model import score

OUT = RESULTS / "literature"


def load_local():
    """accession -> sequence, and gene/protein-name token -> accession."""
    seqs: dict[str, str] = {}
    gene2acc: dict[str, set[str]] = {}
    for name in ("uniprot_lipid_reviewed.tsv", "uniprot_negatives.tsv"):
        df = pd.read_csv(DATA_RAW / name, sep="\t", dtype=str).fillna("")
        for acc, seq, genes in zip(df["Entry"], df["Sequence"], df["Gene Names"]):
            seqs[acc] = seq.strip().upper()
            for g in re.split(r"[ ;]+", genes):
                if len(g) >= 2:
                    gene2acc.setdefault(g.upper(), set()).add(acc)
    mapping = {g: next(iter(v)) for g, v in gene2acc.items() if len(v) == 1}
    print(f"local tables: {len(seqs)} sequences, {len(mapping)} unambiguous gene tokens")
    return seqs, mapping


def resolve(tokens, mapping):
    out = {}
    for t in tokens:
        key = t.upper().replace(" ", "")
        if key in mapping:
            out[t] = mapping[key]
            continue
        pref = [g for g in mapping if g.startswith(key)]
        out[t] = mapping[pref[0]] if len(pref) == 1 else ""
    return out


def main() -> None:
    cand = pd.read_csv(OUT / "candidates.csv")
    seqs, mapping = load_local()
    print("candidates", len(cand), "tokens", cand.token.nunique())
    gmap = resolve(cand.token.unique().tolist(), mapping)
    cand["acc"] = cand.token.map(gmap).fillna("")
    cand = cand[cand.acc != ""]
    print(f"mapped mentions {len(cand)} ({cand.acc.nunique()} proteins, {cand.pmcid.nunique()} papers)")

    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    known = set(zip(sites[sites.label == 1].acc, sites[sites.label == 1].pos))

    rows = []
    for _, r in cand.drop_duplicates(["acc", "pos"]).iterrows():
        seq = seqs.get(r.acc)
        if not seq or r.pos < 1 or r.pos > len(seq):
            continue
        if seq[r.pos - 1] != "C":
            continue
        rows.append({"acc": r.acc, "pos": int(r.pos), "pmcid": r.pmcid, "token": r.token,
                     "window": window(seq, int(r.pos)), "protein_length": len(seq),
                     "in_uniprot": (r.acc, int(r.pos)) in known})
    lit = pd.DataFrame(rows).drop_duplicates(["acc", "pos"])
    print(f"verified cysteine sites: {len(lit)} in {lit.acc.nunique()} proteins "
          f"({lit.in_uniprot.mean():.1%} already in UniProt)")
    lit.to_csv(OUT / "literature_sites.csv", index=False)
    if lit.empty:
        return

    neg_rows = []
    for acc in lit.acc.unique():
        seq = seqs[acc]
        pos_set = set(lit[lit.acc == acc]["pos"])
        cands = [i + 1 for i, a in enumerate(seq) if a == "C" and (i + 1) not in pos_set
                 and (acc, i + 1) not in known]
        if len(cands) > 6:
            step = len(cands) / 6
            cands = [cands[int(i * step)] for i in range(6)]
        for p in cands:
            neg_rows.append({"acc": acc, "pos": p, "pmcid": "", "token": "",
                             "window": window(seq, p), "protein_length": len(seq),
                             "in_uniprot": False})
    neg = pd.DataFrame(neg_rows)
    print(f"same-protein negatives: {len(neg)}")

    df = pd.concat([lit.assign(label=1), neg.assign(label=0)], ignore_index=True)
    df["compartment"] = "unknown"
    pred = score(df[["acc", "pos", "window", "compartment", "protein_length"]])
    df["p"] = pred.p.values
    df.to_csv(OUT / "literature_predictions.csv", index=False)

    thr = json.loads((MODELS / "final_meta.json").read_text(encoding="utf-8"))["threshold_balanced"]
    out = {"n_literature_sites": int(len(lit)), "n_literature_proteins": int(lit.acc.nunique()),
           "n_papers": int(lit.pmcid.nunique()),
           "fraction_in_uniprot": float(lit.in_uniprot.mean()), "threshold": float(thr)}

    nov = df[(df.label == 1) & ~df.in_uniprot]
    negp = df[df.label == 0]
    out["n_novel"] = int(len(nov))
    out["n_novel_proteins"] = int(nov.acc.nunique())
    if len(nov) and len(negp):
        y = np.r_[np.ones(len(nov)), np.zeros(len(negp))]
        p = np.r_[nov.p.values, negp.p.values]
        out["novel_vs_neg_auc"] = float(roc_auc_score(y, p))
        out["novel_vs_neg_aupr"] = float(average_precision_score(y, p))
        out["novel_recall_at_thr"] = float((nov.p.values >= thr).mean())

    # within-protein ranking, controlling for protein-level score offsets
    aucs = []
    for acc, g in df.groupby("acc"):
        pos = g[(g.label == 1) & ~g.in_uniprot]
        negp2 = g[g.label == 0]
        if len(pos) and len(negp2):
            y = np.r_[np.ones(len(pos)), np.zeros(len(negp2))]
            p = np.r_[pos.p.values, negp2.p.values]
            aucs.append((len(pos), roc_auc_score(y, p)))
    if aucs:
        w = np.array([a[0] for a in aucs], dtype=float)
        vals = np.array([a[1] for a in aucs])
        out["within_protein_auc_weighted"] = float((w * vals).sum() / w.sum())
        out["within_protein_auc_mean"] = float(vals.mean())
        out["n_within_protein_tests"] = len(aucs)

    knownp = df[(df.label == 1) & df.in_uniprot]
    if len(knownp):
        out["known_recall_at_thr"] = float((knownp.p.values >= thr).mean())

    # motif sanity check: a PWM from UniProt should enrich literature positives
    tr = sites[sites.window.str.len().eq(WINDOW)]
    fp = np.full((WINDOW, 20), 0.5)
    fn = np.full((WINDOW, 20), 0.5)
    idx = {a: i for i, a in enumerate(AA)}
    for w in tr[tr.label == 1].window:
        for j, c in enumerate(w):
            if c in idx:
                fp[j, idx[c]] += 1
    for w in tr[(tr.label == 0) & (tr.source == "background_cys")].sample(20000, random_state=0).window:
        for j, c in enumerate(w):
            if c in idx:
                fn[j, idx[c]] += 1
    fp /= fp.sum(1, keepdims=True)
    fn /= fn.sum(1, keepdims=True)
    lo = np.log2(fp / fn)

    def pwm(w):
        return float(sum(lo[j, idx[c]] for j, c in enumerate(w) if c in idx))

    df["pwm"] = df.window.map(pwm)
    nov = df[(df.label == 1) & ~df.in_uniprot]
    negp = df[df.label == 0]
    if len(nov) and len(negp):
        y = np.r_[np.ones(len(nov)), np.zeros(len(negp))]
        out["pwm_auc_novel"] = float(roc_auc_score(y, np.r_[nov.pwm.values, negp.pwm.values]))
        out["model_auc_novel"] = float(roc_auc_score(y, np.r_[nov.p.values, negp.p.values]))
    df.to_csv(OUT / "literature_predictions.csv", index=False)
    (OUT / "metrics.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
