"""Leakage-controlled data splits.

Protocol A - `splits_cv.csv`: 5-fold cross-validation over protein clusters.
  Proteins are clustered by 4-mer cosine similarity (default 0.30, an
  approximation of the 30 % identity cut-off used by CD-HIT). Clusters sharing
  an identical 31-aa window are merged, so neither a near-identical protein nor
  a duplicated peptide can appear in both training and test folds.

Protocol B - `split_species.csv`: independent cross-species test.
  Test  = human sites.
  Train = non-human sites after removing (i) every cluster containing a human
          test protein and (ii) every annotation propagated from a human entry
          (evidence "ECO:0000250|UniProtKB:<human acc>"), which would otherwise
          leak the human label into the training set.
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
    chunk = 400
    for start in range(0, len(accs), chunk):
        stop = min(start + chunk, len(accs))
        S = cosine_similarity(X[start:stop], X)
        for i in range(S.shape[0]):
            nz = np.nonzero(S[i, start + i + 1:] >= thr)[0]
            for j in nz:
                uf.union(start + i, start + i + 1 + int(j))
    return uf


def main() -> None:
    sites = pd.read_csv(PROC / "sites.csv", dtype={"organism": str})
    raw = pd.read_csv(RAW / "uniprot_lipid_reviewed.tsv", sep="\t", dtype=str).fillna("")
    neg = pd.read_csv(RAW / "uniprot_negatives.tsv", sep="\t", dtype=str).fillna("")
    seqs = dict(zip(raw["Entry"], raw["Sequence"]))
    seqs.update(dict(zip(neg["Entry"], neg["Sequence"])))
    print(f"sites={len(sites)} unique proteins={sites.acc.nunique()} seqs={len(seqs)}")

    # ---- A. cluster all proteins that carry a positive site ---------------
    pos_accs = sorted(sites.loc[sites.label == 1, "acc"].unique())
    uf = cluster_proteins(pos_accs, seqs)
    clus_of = {a: f"c{uf.find(i):05d}" for i, a in enumerate(pos_accs)}
    print(f"clusters at 4-mer cosine>={SIM_THRESHOLD}: {len(set(clus_of.values()))}")

    parent = {c: c for c in set(clus_of.values())}

    def find(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    win2clus: dict[str, str] = {}
    for _, r in sites[sites.label == 1].iterrows():
        c = clus_of[r.acc]
        if r.window in win2clus:
            a, b = find(c), find(win2clus[r.window])
            if a != b:
                parent[b] = a
        else:
            win2clus[r.window] = c
    clus_of = {a: find(c) for a, c in clus_of.items()}

    sites["cluster"] = sites["acc"].map(clus_of)
    mask_bg = sites["cluster"].isna()
    sites.loc[mask_bg, "cluster"] = "bg_" + sites.loc[mask_bg, "acc"]
    print(f"clusters after window merge: {sites.cluster.nunique()}")

    # ---- balanced folds over clusters -------------------------------------
    stats = (sites.assign(is_pos=(sites.label == 1).astype(int))
             .groupby("cluster").agg(n_pos=("is_pos", "sum")).sort_values("n_pos", ascending=False))
    load = np.zeros(5)
    fold_of = {}
    for c, r in stats.iterrows():
        f = int(np.argmin(load))
        fold_of[c] = f
        load[f] += max(r.n_pos, 1)
    sites["fold"] = sites.cluster.map(fold_of).astype(int)
    sites[["acc", "pos", "label", "cluster", "fold", "source", "organism", "is_human"]].to_csv(
        PROC / "splits_cv.csv", index=False)
    print("positives per fold:", sites[sites.label == 1].groupby("fold").size().to_dict())

    # ---- B. independent cross-species split -------------------------------
    human_prot = set(raw.loc[raw["Organism (ID)"].astype(str).str.startswith("9606"), "Entry"])
    human_prot |= set(neg["Entry"])
    human_test_acc = set(sites.loc[(sites.is_human) & (sites.label == 1), "acc"])
    human_clusters = set(sites.loc[sites.acc.isin(human_test_acc), "cluster"])

    def refs_human(ev: str) -> bool:
        return any(a in human_prot for a in re.findall(r"UniProtKB:([A-Z0-9]+)", str(ev)))

    train_mask = (~sites.is_human) & (~sites.cluster.isin(human_clusters))
    leaky = sites["evidence"].map(refs_human)
    sites["split_species"] = "drop"
    sites.loc[sites.is_human, "split_species"] = "test"
    sites.loc[train_mask & ~leaky, "split_species"] = "train"
    print("species split:", sites.split_species.value_counts().to_dict())
    print("  human positives (test):", int(((sites.split_species == "test") & (sites.label == 1)).sum()))
    print("  train positives:", int(((sites.split_species == "train") & (sites.label == 1)).sum()))
    print("  human-propagated rows excluded:", int((train_mask & leaky).sum()))
    sites[["acc", "pos", "label", "cluster", "split_species", "source", "organism", "is_human"]].to_csv(
        PROC / "split_species.csv", index=False)


if __name__ == "__main__":
    main()
