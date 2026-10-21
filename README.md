# PalmDL

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23455873.svg)](https://doi.org/10.5281/zenodo.23455873)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)

PalmDL predicts S-palmitoylation sites from protein sequence and ships a
leakage-controlled, auditable benchmark. It combines frozen ESM-2 embeddings with a
motif/physicochemical random forest, and provides:

- **cluster-held-out and species-held-out folds** (the canonical benchmark),
- a **nested ensemble** whose members are chosen inside each outer fold,
- a **proteome-wide candidate list** scored with the forest view and empirical q-values,
- a **six-species transfer matrix** (heterogeneous transfer; bacterial N-palmitoyl sites excluded),
- a **sequence-only prediction service** (HTTP + batch, standard library only).

Headline numbers: nested cluster-held-out AUC **0.938** (cluster-bootstrap 95% CI
0.917-0.954), MCC 0.590; human test sites AUC **0.910** (0.890-0.929); strict species
holdout AUC 0.885. Published human AUC 0.900 is context, not a paired comparison.

## Quick start

```bash
pip install -r requirements.txt
python server/app.py --batch server/example.fa --out example.csv --threshold 0.36
python server/app.py            # http://127.0.0.1:8765
```

The service needs model artifacts in `models/` (not versioned). Regenerate them with
`python src/final_model.py` (about 3 min).

## Repository map

| Path | Contents |
|---|---|
| `src/` | pipeline scripts (data, splits, nested ensembles, scan, literature, structure) |
| `data/processed/` | labelled sites and **the canonical fold assignments** |
| `results/` | summary tables, per-fold metrics, shortlists, statistics |
| `figures/`, `manuscript/` | main text and supplementary material |
| `server/`, `deploy/` | sequence-only prediction service and hosting recipes |
| `docs/REPRODUCE.md` | ordered commands, tag v2.1 and file hashes |

## Data and archive

The corrected v2.1 code and processed data live in this repository (tag v2.1). The
earlier v1.0.0 Zenodo record predates the label correction and must not be used for the
current results; a corrected Zenodo version accompanies the manuscript.

## Citation

See `CITATION.cff`; please cite the paper and the archived release.

## Licence

Code MIT (`LICENSE`), derived data CC BY 4.0 (`LICENSE-DATA`).
