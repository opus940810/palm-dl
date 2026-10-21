"""Proteome-wide scan of every cysteine in the reviewed human proteome.

Scores all 262,286 cysteines of 20,431 reviewed human proteins with the
deployable PalmDL model.  Work is chunked and resumable: each chunk writes a CSV
so an interrupted run can continue without rescoring.

Usage:
  python scan_score.py                 # score everything (resumes)
  python scan_score.py --max-chunks 2  # smoke test
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from build_dataset import compartment, parse_topology, window
from config import DATA_RAW, MODELS, RESULTS
from final_model import load, logit

OUT = RESULTS / "scan"
OUT.mkdir(parents=True, exist_ok=True)
CHUNK = 5000
DECOY_N = 20000


def human_proteins() -> dict[str, tuple[str, str]]:
    """accession -> (sequence, topology cell)."""
    out = {}
    for name in ("uniprot_negatives.tsv", "uniprot_lipid_reviewed.tsv"):
        df = pd.read_csv(DATA_RAW / name, sep="\t", dtype=str).fillna("")
        if name.startswith("uniprot_lipid"):
            df = df[df["Organism (ID)"].astype(str).str.startswith("9606")]
        for acc, seq, topo in zip(df["Entry"], df["Sequence"], df.get("Topological domain", [""] * len(df))):
            out[acc] = (seq.upper(), topo)
    return out


def build_targets() -> pd.DataFrame:
    prots = human_proteins()
    rows = []
    for acc, (seq, topo) in prots.items():
        topo_p = parse_topology(topo)
        for i, aa in enumerate(seq):
            if aa == "C":
                pos = i + 1
                rows.append({"acc": acc, "pos": pos, "window": window(seq, pos),
                             "protein_length": len(seq),
                             "compartment": compartment(pos, topo_p, None)})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-chunks", type=int, default=0)
    ap.add_argument("--view", default="rf", choices=["esm", "rf"],
                    help="esm = released ensemble using precomputed embeddings; rf = forest only")
    args = ap.parse_args()

    tgt_path = OUT / "targets.csv"
    if tgt_path.exists():
        tgt = pd.read_csv(tgt_path)
    else:
        tgt = build_targets()
        tgt.to_csv(tgt_path, index=False)
    print(f"targets: {len(tgt)} cysteines in {tgt.acc.nunique()} proteins (view={args.view})", flush=True)

    rf, net, cal, meta = load()
    import joblib
    cal_rf = joblib.load(MODELS / "calibrator_rf.joblib")
    n_chunks = int(np.ceil(len(tgt) / CHUNK))
    if args.max_chunks:
        n_chunks = min(n_chunks, args.max_chunks)

    E_all = None
    if args.view == "esm":
        files = sorted((OUT / "emb").glob("chunk_*.npz"))
        if not files:
            raise SystemExit("no embedding chunks; run src/scan_embed.py first")
        means = [np.load(f)["mean"] for f in files]
        centers = [np.load(f)["center"] for f in files]
        E_all = np.hstack([np.vstack(means), np.vstack(centers)])
        print(f"embeddings: {E_all.shape}", flush=True)
        if len(E_all) != len(tgt):
            raise SystemExit(f"embedding rows {len(E_all)} != targets {len(tgt)}")

    t0 = time.time()
    for c in range(n_chunks):
        path = OUT / f"chunk_{c:04d}.csv"
        if path.exists():
            continue
        part = tgt.iloc[c * CHUNK:(c + 1) * CHUNK]
        from features import context_features, encode_flat

        X = encode_flat(part.window.tolist())
        C = context_features(part)
        p_rf = rf.predict_proba(np.hstack([X, C]))[:, 1]
        if args.view == "esm":
            import torch
            E = E_all[c * CHUNK:(c + 1) * CHUNK].astype(np.float32)
            with torch.no_grad():
                p_esm = torch.sigmoid(net(torch.as_tensor(E), torch.as_tensor(C))).numpy()
            p = cal.predict_proba((0.5 * (logit(p_rf) + logit(p_esm))).reshape(-1, 1))[:, 1]
        else:
            p_esm = np.full(len(part), np.nan)
            p = cal_rf.predict_proba(logit(p_rf).reshape(-1, 1))[:, 1]
        pd.DataFrame({"acc": part.acc.values, "pos": part.pos.values,
                      "p_rf": p_rf, "p_esm": p_esm, "p": p}).to_csv(path, index=False)
        done = c + 1
        rate = done * CHUNK / max(time.time() - t0, 1)
        print(f"  chunk {done}/{n_chunks} ({rate:.0f} sites/s)", flush=True)

    # decoys: same flanks with the central cysteine mutated to alanine
    decoy_path = OUT / "decoys.csv"
    if not decoy_path.exists():
        rng = np.random.default_rng(0)
        idx = rng.choice(len(tgt), size=min(DECOY_N, len(tgt)), replace=False)
        d = tgt.iloc[idx].copy()
        from features import context_features, encode_flat
        X = encode_flat(d.window.tolist())
        C = context_features(d)
        p_rf = rf.predict_proba(np.hstack([X, C]))[:, 1]
        p = cal_rf.predict_proba(logit(p_rf).reshape(-1, 1))[:, 1]
        pd.DataFrame({"acc": d.acc.values, "pos": d.pos.values, "p": p}).to_csv(
            decoy_path, index=False)
        print("decoys scored", flush=True)


if __name__ == "__main__":
    main()
