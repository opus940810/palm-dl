"""Extract ESM-2 embeddings for every 31-aa window in the dataset.

Outputs (float16, rows aligned with data/processed/sites.csv order):
  data/embeddings/esm2_mean.npy    mean-pooled representation (N, 480)
  data/embeddings/esm2_center.npy  representation at the central cysteine (N, 480)
"""
from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import torch

import esm
from config import DATA_PROC, EMB, ROOT, WINDOW

torch.set_num_threads(12)
CKPT = ROOT.parent / ".cache" / "esm2_t12_35M_UR50D.pt"


def load_esm(path):
    """Load a fair-esm checkpoint (pickled objects) under torch>=2.6 defaults."""
    orig = torch.load

    def _load(*a, **k):
        k["weights_only"] = False  # checkpoint is the official Facebook weights file
        return orig(*a, **k)

    torch.load = _load
    try:
        return esm.pretrained.load_model_and_alphabet_local(str(path))
    finally:
        torch.load = orig


def main() -> None:
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    sites = sites[sites.window.str.len() == WINDOW].reset_index(drop=True)
    windows = sites.window.tolist()
    uniq = list(dict.fromkeys(windows))  # preserve order, deduplicate
    print(f"{len(windows)} rows, {len(uniq)} unique windows", flush=True)

    model, alphabet = load_esm(CKPT)
    model.eval()
    bc = alphabet.get_batch_converter()
    dim = model.embed_dim

    mean = np.zeros((len(uniq), dim), dtype=np.float16)
    center = np.zeros((len(uniq), dim), dtype=np.float16)
    bs = 64
    t0 = time.time()
    with torch.no_grad():
        for start in range(0, len(uniq), bs):
            chunk = uniq[start:start + bs]
            _, _, toks = bc([(f"s{i}", s) for i, s in enumerate(chunk)])
            out = model(toks, repr_layers=[model.num_layers], return_contacts=False)
            rep = out["representations"][model.num_layers]  # (B, L+2, dim)
            mean[start:start + len(chunk)] = rep[:, 1:-1].mean(dim=1).numpy().astype(np.float16)
            center[start:start + len(chunk)] = rep[:, 1 + WINDOW // 2].numpy().astype(np.float16)
            if start % (bs * 20) == 0:
                done = start + len(chunk)
                print(f"  {done}/{len(uniq)} ({time.time() - t0:.0f}s)", flush=True)

    idx = {w: i for i, w in enumerate(uniq)}
    order = np.array([idx[w] for w in windows])
    np.save(EMB / "esm2_mean.npy", mean[order])
    np.save(EMB / "esm2_center.npy", center[order])
    print(f"saved {mean[order].shape} in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
