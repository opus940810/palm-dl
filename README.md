# PalmDL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23455873.svg)](https://doi.org/10.5281/zenodo.23455873)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

PalmDL predicts S-palmitoylation sites from protein sequence and ships a
leakage-controlled benchmark. It combines frozen ESM-2 embeddings with a
motif/physicochemical random forest, and provides:

- **cluster-held-out and species-held-out folds** (the canonical benchmark),
- a **proteome-wide, FDR-controlled candidate list** for the human proteome,
- a **cross-species transfer matrix** (eukaryotic code shared; bacterial lipidation distinct),
- a **structure-aware variant** using AlphaFold pLDDT/packing (graceful fallback),
- a **prediction service** (HTTP + batch, standard library only).

Headline numbers: cluster-held-out AUC **0.956** (95% CI 0.953–0.959); human test
sites AUC **0.901** (0.885–0.915), statistically indistinguishable from the best
published human value 0.900 (bootstrap P = 0.46); strict species holdout AUC 0.822.

## Quick start

```bash
pip install -r requirements.txt
python server/app.py --batch server/example.fa --out example.csv --threshold 0.5
python server/app.py            # http://127.0.0.1:8765
```

The service needs model artifacts in `models/` (not versioned). Regenerate them
with `python src/final_model.py` (about 3 min) and, for the structure path,
`python src/structure_hybrid.py` (about 10 min).

## Repository map

| Path | Contents |
|---|---|
| `src/` | 25 pipeline scripts (data, splits, models, scan, literature, structure) |
| `data/processed/` | labelled sites and **the canonical fold assignments** |
| `results/` | summary tables, per-fold metrics, shortlists, statistics |
| `figures/`, `manuscript/` | main text and supplementary material |
| `server/`, `deploy/` | prediction service and hosting recipes |
| `docs/REPRODUCE.md` | ordered commands reproducing every number |

## Data and archive

This repository is archived at Zenodo (DOI 10.5281/zenodo.23455873). The exact
UniProt 2026_03 snapshot, out-of-fold predictions, the proteome-wide score table
and the complete release archive are attached to the GitHub release `v1.0.0` as
`palmdl-v1.0.0.zip`, because they are too large for git history. Raw data are
CC BY 4.0 (UniProt), AlphaFold models CC BY 4.0.

## Citation and priority

See `CITATION.cff`; please cite the paper and the Zenodo DOI. The release timeline
(GitHub commit 2026-10-18, Zenodo DOI 10.5281/zenodo.23455873) is recorded in
[`PRIORITY.md`](PRIORITY.md); if you evaluate a new method, please use the canonical
folds in `data/processed/splits_cv.csv` so that numbers stay comparable.

## Licence

Code MIT (`LICENSE`), derived data CC BY 4.0 (`LICENSE-DATA`).
