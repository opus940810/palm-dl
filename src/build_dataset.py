"""Build the site-level dataset from the raw UniProt TSV exports.

Label definition
----------------
positive : cysteine annotated with a 'Lipidation' feature whose note is
           "S-palmitoyl cysteine".
negative : (i)  other cysteines of the same protein that are not annotated and
                whose +/-15 window does not overlap an annotated site
                ("hard" / same-protein negatives);
           (ii) cysteines of reviewed human proteins that carry no lipidation
                feature at all ("background" negatives).

Every record keeps the evidence code so that models can be trained on
high-confidence sites only and evaluated on experimentally supported ones.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
PROC.mkdir(parents=True, exist_ok=True)

HALF = 15  # window = 31 aa, same as CSS-Palm 4.0 / Deep-Palm

LIPID_RE = re.compile(r'LIPID (\d+)(?:\.\.(\d+))?;\s*/note="([^"]*)";(?:\s*/evidence="([^"]*)")?')
TOPO_RE = re.compile(r'TOPO_DOM (\d+)(?:\.\.(\d+))?;\s*/note="([^"]*)";(?:\s*/evidence="([^"]*)")?')
SIGNAL_RE = re.compile(r'SIGNAL (\d+)\.\.(\d+)')

EXPERIMENTAL = ("ECO:0000269", "ECO:0007744", "ECO:0007005", "ECO:0006056")


def parse_lipidation(cell: str) -> list[tuple[int, str, str]]:
    """Return [(position, note, evidence), ...] for one TSV cell."""
    out = []
    if not isinstance(cell, str):
        return out
    for m in LIPID_RE.finditer(cell):
        start = int(m.group(1))
        note = (m.group(3) or "").strip()
        ev = (m.group(4) or "").strip()
        out.append((start, note, ev))
    return out


def parse_topology(cell: str) -> list[tuple[int, int, str]]:
    out = []
    if not isinstance(cell, str):
        return out
    for m in TOPO_RE.finditer(cell):
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        note = (m.group(3) or "").strip().lower()
        out.append((start, end, note))
    return out


def compartment(pos: int, topo: list[tuple[int, int, str]], signal: tuple[int, int] | None) -> str:
    if signal and signal[0] <= pos <= signal[1]:
        return "signal_peptide"
    for start, end, note in topo:
        if start <= pos <= end:
            if "cytoplasm" in note:
                return "cytoplasmic"
            if "extracellular" in note:
                return "extracellular"
            if "lumen" in note:
                return "lumenal"
            return "other"
    return "unknown"


def window(seq: str, pos: int, half: int = HALF) -> str:
    """1-based pos; zero-padded to fixed length 2*half+1."""
    s = pos - 1
    lo, hi = s - half, s + half + 1
    left_pad, right_pad = max(0, -lo), max(0, hi - len(seq))
    core = seq[max(0, lo): min(len(seq), hi)]
    return "X" * left_pad + core + "X" * right_pad


def is_experimental(ev: str) -> bool:
    return any(code.split("|")[0] == c for code in ev.split(",") for c in EXPERIMENTAL)


def build() -> None:
    lip = pd.read_csv(RAW / "uniprot_lipid_reviewed.tsv", sep="\t", dtype=str).fillna("")
    neg = pd.read_csv(RAW / "uniprot_negatives.tsv", sep="\t", dtype=str).fillna("")
    print(f"lipid entries: {len(lip)}, background entries: {len(neg)}")

    rows = []
    prot_rows = []
    for _, r in lip.iterrows():
        seq = r["Sequence"].strip().upper()
        acc = r["Entry"]
        feats = parse_lipidation(r.get("Lipidation", ""))
        sites = [(p, note, ev) for p, note, ev in feats if "palmitoyl" in note.lower()]
        if not sites:
            continue
        seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", seq)
        topo = parse_topology(r.get("Topological domain", ""))
        sig = None
        m = SIGNAL_RE.search(r.get("Signal peptide", "") or "")
        if m:
            sig = (int(m.group(1)), int(m.group(2)))
        site_pos = {p for p, _, _ in sites}
        prot_rows.append(
            dict(acc=acc, organism=r.get("Organism (ID)", ""), protein=r.get("Protein names", "")[:120],
                 gene=r.get("Gene Names", ""), length=len(seq), n_sites=len(site_pos))
        )
        for p, note, ev in sites:
            if p < 1 or p > len(seq):
                continue
            rows.append(dict(acc=acc, pos=p, label=1, residue=seq[p - 1],
                             source="uniprot_exp" if is_experimental(ev) else "uniprot_inferred",
                             evidence=ev, compartment=compartment(p, topo, sig),
                             window=window(seq, p), protein_length=len(seq),
                             organism=r.get("Organism (ID)", ""), gene=r.get("Gene Names", "")))
        # hard negatives: unannotated cysteines of the same protein
        cand = [i + 1 for i, aa in enumerate(seq) if aa == "C" and (i + 1) not in site_pos]
        cand = [p for p in cand if all(abs(p - s) > HALF for s in site_pos)]
        # deterministic, spread-out subset of at most 4 per protein
        if len(cand) > 4:
            step = len(cand) / 4
            cand = [cand[int(i * step)] for i in range(4)]
        for p in cand:
            rows.append(dict(acc=acc, pos=p, label=0, residue="C", source="same_protein_cys",
                             evidence="", compartment=compartment(p, topo, sig),
                             window=window(seq, p), protein_length=len(seq),
                             organism=r.get("Organism (ID)", ""), gene=r.get("Gene Names", "")))

    # background negatives (human proteins without any lipidation annotation)
    n_bg = 0
    for _, r in neg.iterrows():
        seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", r["Sequence"].strip().upper())
        topo = parse_topology(r.get("Topological domain", ""))
        sig_m = SIGNAL_RE.search(r.get("Signal peptide", "") or "")
        sig = (int(sig_m.group(1)), int(sig_m.group(2))) if sig_m else None
        cand = [i + 1 for i, aa in enumerate(seq) if aa == "C"]
        if not cand:
            continue
        if len(cand) > 3:
            step = len(cand) / 3
            cand = [cand[int(i * step)] for i in range(3)]
        for p in cand:
            rows.append(dict(acc=r["Entry"], pos=p, label=0, residue="C", source="background_cys",
                             evidence="", compartment=compartment(p, topo, sig),
                             window=window(seq, p), protein_length=len(seq),
                             organism=r.get("Organism (ID)", ""), gene=r.get("Gene Names", "")))
            n_bg += 1

    df = pd.DataFrame(rows)
    df = df.drop_duplicates(subset=["acc", "pos"]).reset_index(drop=True)
    df["is_human"] = df["organism"].astype(str).str.startswith("9606")
    df.to_csv(PROC / "sites.csv", index=False)
    pd.DataFrame(prot_rows).drop_duplicates("acc").to_csv(PROC / "proteins.csv", index=False)

    print(df.groupby(["source"]).size())
    print("\nlabel x source\n", pd.crosstab(df["label"], df["source"]))
    print("\npositives by organism (top 10)\n", df[df.label == 1]["organism"].value_counts().head(10))
    print("\ncompartment of positives\n", df[df.label == 1]["compartment"].value_counts())
    print("\nunique human positive sites:", df[(df.label == 1) & df.is_human].shape[0])
    print("unique experimental positive sites:", df[(df.label == 1) & (df.source == "uniprot_exp")].shape[0])
    print("background negatives:", n_bg, "same-protein negatives:", (df.source == "same_protein_cys").sum())


if __name__ == "__main__":
    build()
