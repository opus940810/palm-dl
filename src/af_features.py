"""AlphaFold structural features for human palmitoylation sites.

Downloads AlphaFold DB models (via the prediction API, which carries the current
model URL) for human proteins and computes per-site features from the PDB:

  plddt_site    pLDDT of the cysteine C-alpha
  plddt_window  mean pLDDT over residues within +/-5 of the site
  packing       number of C-alpha atoms within 10 A of the site C-alpha
                (excluding sequence neighbours |i-j| < 4)

Missing residues are tolerated; a site is scored only if a C-alpha within +/-2
residues of the annotated position exists.
"""
from __future__ import annotations

import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from config import DATA_RAW, RESULTS

OUT = RESULTS / "structure"
CACHE = OUT / "af"
CACHE.mkdir(parents=True, exist_ok=True)
H = {"User-Agent": "Mozilla/5.0"}
MAX_LEN = 700


def human_proteins() -> dict[str, str]:
    seqs = {}
    for name in ("uniprot_negatives.tsv", "uniprot_lipid_reviewed.tsv"):
        df = pd.read_csv(DATA_RAW / name, sep="\t", dtype=str).fillna("")
        if name.startswith("uniprot_lipid"):
            df = df[df["Organism (ID)"].astype(str).str.startswith("9606")]
        for acc, seq in zip(df["Entry"], df["Sequence"]):
            seqs[acc] = seq.upper()
    return seqs


def fetch(acc: str) -> str | None:
    pdb = CACHE / f"{acc}.pdb"
    if pdb.exists() and pdb.stat().st_size > 1000:
        return str(pdb)
    try:
        meta = json.load(urllib.request.urlopen(
            urllib.request.Request(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", headers=H),
            timeout=60))
        if not meta:
            return None
        url = meta[0]["pdbUrl"]
        data = urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=120).read()
        pdb.write_bytes(data)
        return str(pdb)
    except Exception:
        return None


def parse_pdb(path: str) -> dict[int, tuple[float, float, float, float]]:
    """residue number -> (x, y, z, pLDDT) for C-alpha atoms."""
    out = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                resnum = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                b = float(line[60:66])
            except ValueError:
                continue
            out[resnum] = (x, y, z, b)
    return out


def features(sites: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for acc, g in sites.groupby("acc"):
        path = CACHE / f"{acc}.pdb"
        if not path.exists():
            continue
        ca = parse_pdb(str(path))
        if not ca:
            continue
        nums = np.array(sorted(ca))
        xyz = np.array([ca[n][:3] for n in nums])
        pld = np.array([ca[n][3] for n in nums])
        order = {n: i for i, n in enumerate(nums)}
        for _, s in g.iterrows():
            pos = int(s.pos)
            # nearest modelled residue within +/-2
            cand = [n for n in (pos, pos - 1, pos + 1, pos - 2, pos + 2) if n in order]
            if not cand:
                continue
            p = cand[0]
            i = order[p]
            x = xyz[i]
            d = np.linalg.norm(xyz - x, axis=1)
            neigh = (d < 10.0) & (np.abs(np.arange(len(nums)) - i) >= 4)
            win = (np.abs(nums - p) <= 5)
            rows.append({"acc": acc, "pos": pos, "modelled": p,
                         "plddt_site": float(pld[i]),
                         "plddt_window": float(pld[win].mean()),
                         "packing": int(neigh.sum()),
                         "n_residues": len(nums)})
    return pd.DataFrame(rows)


def main() -> None:
    seqs = human_proteins()
    sites = pd.read_csv(DATA_PROC_SITES(), dtype={"organism": str}, low_memory=False)
    sites = sites[sites.is_human & sites.window.str.len().eq(31)]
    pos_acc = set(sites[sites.label == 1].acc)
    bg_acc = [a for a in seqs if a not in pos_acc and len(seqs[a]) <= MAX_LEN]
    rng = np.random.default_rng(0)
    bg_acc = list(rng.choice(bg_acc, size=min(1200, len(bg_acc)), replace=False))
    wanted = [a for a in list(pos_acc) + bg_acc if a in seqs and len(seqs[a]) <= MAX_LEN]
    print(f"human proteins to fetch: {len(wanted)} ({len(pos_acc)} with sites)", flush=True)

    with ThreadPoolExecutor(max_workers=6) as ex:
        got = list(ex.map(fetch, wanted))
    ok = sum(1 for g in got if g)
    print(f"models available: {ok}/{len(wanted)}", flush=True)

    sub = sites[sites.acc.isin({a for a, g in zip(wanted, got) if g})]
    feats = features(sub)
    feats.to_csv(OUT / "site_structure_features.csv", index=False)
    print(f"structural features for {len(feats)} sites in {feats.acc.nunique()} proteins")


def DATA_PROC_SITES():
    from config import DATA_PROC
    return DATA_PROC / "sites.csv"


if __name__ == "__main__":
    main()
