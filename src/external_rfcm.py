"""External validation on the RFCM-PALM mouse synaptic hold-out data.

The RFCM-PALM repository (github.com/anupgth/RFCM-PALM) distributes the
PANIMoni-derived mouse synaptic S-palmitoylation data as 39-aa windows with a
central cysteine, split into training and hold-out sets.  We use the hold-out
sets only, convert the 39-mer to our 31-mer, and score with the deployable
PalmDL model.  Sites whose (accession, position) already occur in our UniProt
training table are reported separately ("seen" vs "unseen").
"""
from __future__ import annotations

import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from config import DATA_PROC, DATA_RAW, MODELS, RESULTS, WINDOW
from final_model import score
from metrics import metrics_at

EXT = Path("C:/Users/Tyrese/Documents/deepseek-harness/default-workspace/.cache/rfcm/ext")
OUT = RESULTS / "external"
OUT.mkdir(exist_ok=True)
HALF39 = 19
OFF = HALF39 - WINDOW // 2  # 4


def parse(path: Path, label: int) -> pd.DataFrame:
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        acc, pos, win = parts[0], int(parts[1]), parts[2]
        if len(win) != 2 * HALF39 + 1 or win[HALF39] != "C":
            continue
        rows.append({"acc": acc, "pos": pos, "window": win[OFF:OFF + WINDOW], "label": label,
                     "source": "rfcm_ho"})
    return pd.DataFrame(rows)


def fetch_lengths(accs: list[str]) -> dict[str, int]:
    out = {}
    for i in range(0, len(accs), 50):
        chunk = accs[i:i + 50]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": q, "format": "tsv", "size": 500, "fields": "accession,length"})
        try:
            txt = urllib.request.urlopen(urllib.request.Request(
                url, headers={"User-Agent": "Mozilla/5.0"}), timeout=180).read().decode()
            for line in txt.splitlines()[1:]:
                f = line.split("\t")
                if len(f) >= 2:
                    out[f[0]] = int(f[1])
        except Exception as exc:  # noqa: BLE001
            print("length fetch failed", chunk[:2], repr(exc)[:80])
    return out


def main() -> None:
    frames = []
    for subset in ("Combined", "Male", "Female"):
        p = parse(EXT / "TR-Data" / subset / "HO-Test" / "TestSeq_pos.txt", 1)
        n = parse(EXT / "TR-Data" / subset / "HO-Test" / "TestSeq_neg.txt", 0)
        d = pd.concat([p, n], ignore_index=True)
        d["subset"] = subset
        frames.append(d)
        print(subset, len(p), "pos", len(n), "neg")
    ext = pd.concat(frames, ignore_index=True)

    # protein length for context features: local raw tables first, then UniProt
    seqs = {}
    for f in ("uniprot_lipid_reviewed.tsv", "uniprot_negatives.tsv"):
        raw = pd.read_csv(DATA_RAW / f, sep="\t", dtype=str).fillna("")
        seqs.update(dict(zip(raw["Entry"], raw["Sequence"])))
    lengths = {a: len(s) for a, s in seqs.items()}
    missing = [a for a in ext.acc.unique() if a not in lengths]
    if missing:
        print("fetching lengths for", len(missing), "accessions")
        lengths.update(fetch_lengths(missing))
    ext["protein_length"] = ext.acc.map(lengths).fillna(ext.pos + 20).astype(int)
    ext["compartment"] = "unknown"

    # seen / unseen split against our UniProt table
    sites = pd.read_csv(DATA_PROC / "sites.csv", dtype={"organism": str}, low_memory=False)
    key = set(zip(sites.acc, sites.pos))
    ext["seen"] = [bool((a, p) in key) for a, p in zip(ext.acc, ext.pos)]
    print("seen fraction:", round(ext.seen.mean(), 3), "| positives seen:",
          round(ext[ext.label == 1].seen.mean(), 3))

    pred = score(ext[["acc", "pos", "window", "compartment", "protein_length"]])
    ext = ext.assign(p_rf=pred.p_rf.values, p_esm=pred.p_esm.values, p=pred.p.values)

    meta = json.loads((MODELS / "final_meta.json").read_text(encoding="utf-8"))
    thr_bal = meta["threshold_balanced"]

    rows = []
    for subset in ("Combined", "Male", "Female"):
        for unseen_only in (False, True):
            d = ext[ext.subset == subset]
            if unseen_only:
                d = d[~d.seen]
            y = d.label.values
            if len(np.unique(y)) < 2 or len(d) < 30:
                continue
            m = metrics_at(y, d.p.values, thr_bal)
            # sensitivity at 95 % specificity, threshold taken from our validation negatives
            rows.append({
                "subset": subset, "unseen_only": unseen_only,
                "n_pos": int(y.sum()), "n_neg": int((y == 0).sum()),
                "auc": roc_auc_score(y, d.p.values),
                "aupr": average_precision_score(y, d.p.values),
                "mcc_at_bal_thr": m["mcc"], "sn_at_bal_thr": m["sn"], "sp_at_bal_thr": m["sp"],
                "threshold": thr_bal,
            })
            print(rows[-1])

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "external_rfcm_metrics.csv", index=False)
    ext.to_csv(OUT / "external_rfcm_predictions.csv", index=False)
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
