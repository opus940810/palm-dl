"""Leakage-controlled data splits (v2, after the 2026-10-19 review).

Protocol A - `splits_cv.csv`: 5-fold cross-validation.
  *Every* protein in the dataset (positives and background negatives) is
  clustered by 4-mer cosine similarity (0.30, a coarse homology proxy that is
  NOT 30 % sequence identity). Clusters that share an identical 31-aa window -
  including windows shared between a positive protein and a background protein -
  are merged, so no exact peptide occurs in two folds. Folds are balanced by
  positive count.

Protocol B - `split_species.csv`: independent cross-species test.
  Test  = human positives and human same-protein negatives.
  Train = non-human rows outside every cluster that contains a human protein,
          after removing annotations propagated from a human entry
          (evidence "ECO:0000250|UniProtKB:<human acc>").
  Background negatives are human by construction and are therefore excluded
  from both sides of the species split (v1 wrongly placed them in the test set,
  which made "no human data in training" only half true).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics.pairwise import cosine_similarity

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"

K = 4
SIM_THRESHOLD = 0.30
N_FOLDS = 5


def build_kmer(seqs: list[str], k: int = K) -> sparse.csr_matrix:
    idx = {b: i for i, b in enumerate("ACDEFGHIKLMNPQRSTVWY")}
    indptr, indices, data = [0], [], []
    for s in seqs:
        v = np.zeros(20 ** k, dtype=np.float32)
        for i in range(len(s) - k + 1):
            j, ok = 0, True
            for ch in s[i:i + k]:
                if ch not in idx:
                    ok = False
                    break
                j = j * 20 + idx[ch]
            if ok:
                v[j] += 1
        v = np.sqrt(v)
        nz = np.nonzero(v)[0]
        indices.extend(nz.tolist())
        data.extend(v[nz].tolist())
        indptr.append(len(indices))
    return sparse.csr_matrix((data, indices, indptr), shape=(len(seqs), 20 ** k))


class UnionFind:
    def __init__(self, n: int):
        self.p = list(range(n))

    def find(self, a: int) -> int:
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def cluster_proteins(accs: list[str], seqs: dict[str, str], thr: float = SIM_THRESHOLD) -> UnionFind:
    X = build_kmer([seqs.get(a, "") for a in accs])
    uf = UnionFind(len(accs))
    chunk = 300
    for start in range(0, len(accs), chunk):
        stop = min(start + chunk, len(accs))
        S = cosine_similarity(X[start:stop], X)
        for i in range(S.shape[0]):
            row = S[i]
            nz = np.nonzero(row >= thr)[0]
            for j in nz:
                if j > start + i:
                    uf.union(start + i, int(j))
    return uf


def main() -> None:
    sites = pd.read_csv(PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    raw = pd.read_csv(RAW / "uniprot_lipid_reviewed.tsv", sep="\t", dtype=str).fillna("")
    neg = pd.read_csv(RAW / "uniprot_negatives.tsv", sep="\t", dtype=str).fillna("")
    seqs = dict(zip(raw["Entry"], raw["Sequence"]))
    seqs.update(dict(zip(neg["Entry"], neg["Sequence"])))
    print(f"sites={len(sites)} unique proteins={sites.acc.nunique()} seqs={len(seqs)}", flush=True)

    # ---- A. cluster every protein that appears in the dataset --------------
    accs = sorted(sites.acc.unique())
    print(f"clustering {len(accs)} proteins at 4-mer cosine >= {SIM_THRESHOLD} ...", flush=True)
    uf = cluster_proteins(accs, seqs)
    group = {a: f"c{uf.find(i):05d}" for i, a in enumerate(accs)}
    print(f"clusters before window merge: {len(set(group.values()))}", flush=True)

    # merge groups that share an identical window (positives and negatives alike)
    parent = {g: g for g in set(group.values())}

    def find(g: str) -> str:
        while parent[g] != g:
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g

    win2grp: dict[str, str] = {}
    merged = 0
    for acc, w in zip(sites.acc, sites.window):
        g = find(group[acc])
        if w in win2grp:
            a, b = g, find(win2grp[w])
            if a != b:
                parent[b] = a
                merged += 1
        else:
            win2grp[w] = g
    group = {a: find(g) for a, g in group.items()}
    sites["cluster"] = sites.acc.map(group)
    print(f"window merges: {merged}; final clusters: {sites.cluster.nunique()}", flush=True)

    # ---- balanced folds over clusters -------------------------------------
    stats = (sites.assign(is_pos=(sites.label == 1).astype(int))
             .groupby("cluster").agg(n_pos=("is_pos", "sum")).sort_values("n_pos", ascending=False))
    load = np.zeros(N_FOLDS)
    fold_of = {}
    for c, r in stats.iterrows():
        f = int(np.argmin(load))
        fold_of[c] = f
        load[f] += max(r.n_pos, 1)
    sites["fold"] = sites.cluster.map(fold_of).astype(int)
    sites[["acc", "pos", "label", "cluster", "fold", "source", "organism", "is_human"]].to_csv(
        PROC / "splits_cv.csv", index=False)
    print("positives per fold:", sites[sites.label == 1].groupby("fold").size().to_dict())

    # integrity: no exact window may appear in two folds
    wf = sites.groupby("window").fold.nunique()
    bad = int((wf > 1).sum())
    print(f"windows crossing folds: {bad} (must be 0)")
    pos = sites[sites.label == 1]
    print(f"positive windows crossing folds: {int((pos.groupby('window').fold.nunique() > 1).sum())}")

    # ---- B. independent cross-species split -------------------------------
    human_prot = set(raw.loc[raw["Organism (ID)"].astype(str).str.startswith("9606"), "Entry"])
    human_prot |= set(neg["Entry"])
    human_clusters = set(sites.loc[sites.is_human, "cluster"])

    def refs_human(ev: str) -> bool:
        return any(a in human_prot for a in re.findall(r"UniProtKB:([A-Z0-9]+)", str(ev)))

    bg = sites.source == "background_cys"
    train_mask = (~sites.is_human) & (~sites.cluster.isin(human_clusters)) & (~bg)
    leaky = sites["evidence"].map(refs_human)
    sites["split_species"] = "drop"
    sites.loc[sites.is_human & (~bg), "split_species"] = "test"
    sites.loc[train_mask & ~leaky, "split_species"] = "train"
    print("species split:", sites.split_species.value_counts().to_dict())
    print("  human positives (test):", int(((sites.split_species == "test") & (sites.label == 1)).sum()))
    print("  train positives:", int(((sites.split_species == "train") & (sites.label == 1)).sum()))
    print("  human-propagated rows excluded:", int((train_mask & leaky).sum()))
    sites[["acc", "pos", "label", "cluster", "split_species", "source", "organism", "is_human"]].to_csv(
        PROC / "split_species.csv", index=False)


if __name__ == "__main__":
    main()
