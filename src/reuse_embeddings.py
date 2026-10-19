"""Reuse the previous ESM embeddings by mapping 31-aa windows, then embed only misses.

The v1 cache (data/embeddings/esm2_*.npy) is aligned with the archived
mixed-label sites table. Windows are pure sequence strings, so identical windows
produce identical embeddings regardless of the label set: we build a
window -> row map from the archive and fill the new table from it. Only windows
that never occurred before are embedded (usually a small minority).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
EMB = ROOT / "data" / "embeddings"
ARCHIVE = PROC / "_archive" / "sites_v1_mixed_labels.csv"
WINDOW = 31


def main() -> None:
    new = pd.read_csv(PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    new = new[new.window.str.len() == WINDOW].reset_index(drop=True)
    old = pd.read_csv(ARCHIVE, low_memory=False)
    old = old[old.window.str.len() == WINDOW].reset_index(drop=True)
    mean_old = np.load(EMB / "esm2_mean.npy", mmap_mode="r")
    center_old = np.load(EMB / "esm2_center.npy", mmap_mode="r")
    assert len(old) == len(mean_old) == len(center_old), (len(old), len(mean_old), len(center_old))

    row_of: dict[str, int] = {}
    for i, w in enumerate(old.window):
        row_of.setdefault(w, i)

    mean = np.zeros((len(new), mean_old.shape[1]), dtype=np.float16)
    center = np.zeros((len(new), center_old.shape[1]), dtype=np.float16)
    missing: list[int] = []
    hit = 0
    for i, w in enumerate(new.window):
        j = row_of.get(w)
        if j is None:
            missing.append(i)
        else:
            mean[i] = mean_old[j]
            center[i] = center_old[j]
            hit += 1
    print(f"rows={len(new)}  reused={hit} ({hit/len(new):.1%})  missing={len(missing)}")

    if missing:
        miss_path = EMB / "missing_windows.npy"
        np.save(miss_path, np.array([new.window[i] for i in missing], dtype=object))
        print(f"saved {len(missing)} windows needing embedding -> {miss_path}")
    np.save(EMB / "esm2_mean_reused.npy", mean)
    np.save(EMB / "esm2_center_reused.npy", center)
    print("saved reused arrays")


if __name__ == "__main__":
    sys.exit(main())
