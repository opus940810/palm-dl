"""Embed the proteome-wide scan windows in resumable chunks.

Writes results/scan/emb/chunk_XXXX.npz with mean/center float16 arrays for the
window slice [start, stop). Re-running skips chunks that already exist, so a
stalled or killed process can be resumed without recomputing anything.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import esm
from config import RESULTS, ROOT

OUT = RESULTS / "scan" / "emb"
OUT.mkdir(parents=True, exist_ok=True)
CKPT = ROOT.parent / ".cache" / "esm2_t12_35M_UR50D.pt"
CHUNK = 10_000
BS = 64


def load_esm(path: Path):
    orig = torch.load

    def _load(*a, **k):
        k["weights_only"] = False
        return orig(*a, **k)

    torch.load = _load
    try:
        return esm.pretrained.load_model_and_alphabet_local(str(path))
    finally:
        torch.load = orig


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--stop", type=int, default=-1)
    args = ap.parse_args()
    torch.set_num_threads(args.threads)

    tgt = pd.read_csv(RESULTS / "scan" / "targets.csv")
    windows = tgt.window.tolist()
    n = len(windows)
    stop_all = n if args.stop < 0 else min(args.stop, n)
    print(f"range {args.start}-{stop_all} of {n} windows, chunk={CHUNK}", flush=True)

    model, alphabet = load_esm(CKPT)
    model.eval()
    bc = alphabet.get_batch_converter()
    dim = model.embed_dim

    for start in range(args.start, stop_all, CHUNK):
        stop = min(start + CHUNK, stop_all)
        path = OUT / f"chunk_{start:06d}.npz"
        if path.exists():
            continue
        t0 = time.time()
        mean = np.zeros((stop - start, dim), dtype=np.float16)
        center = np.zeros((stop - start, dim), dtype=np.float16)
        with torch.no_grad():
            for i in range(start, stop, BS):
                chunk = windows[i:min(i + BS, stop)]
                _, _, toks = bc([(str(j), s) for j, s in enumerate(chunk)])
                out = model(toks, repr_layers=[model.num_layers], return_contacts=False)
                rep = out["representations"][model.num_layers]
                mean[i - start:i - start + len(chunk)] = rep[:, 1:-1].mean(dim=1).numpy().astype(np.float16)
                center[i - start:i - start + len(chunk)] = rep[:, 1 + 31 // 2].numpy().astype(np.float16)
        np.savez_compressed(path, mean=mean, center=center)
        print(f"chunk {start}-{stop} in {time.time() - t0:.0f}s", flush=True)
    print("all chunks done", flush=True)


if __name__ == "__main__":
    main()
