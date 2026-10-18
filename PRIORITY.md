# Priority record

This repository is a timestamped public record of the PalmDL benchmark, model and
analysis. It was released before journal submission so that the priority of the
work is documented and citable, and so that the benchmark folds are available to
anyone who wants to reproduce or extend the results.

## Timeline (UTC)

| Date | Event | Evidence |
|---|---|---|
| 2026-10-18 16:31 | GitHub repository created and first commit pushed | `created_at` / commit `aeeffba` |
| 2026-10-18 16:47 | Release `v1.0.0` published (tag, full archive asset) | GitHub Releases |
| 2026-10-18 | Zenodo record published, DOI `10.5281/zenodo.23455873` | https://doi.org/10.5281/zenodo.23455873 |
| 2026-10-18 16:50 | Metadata commit with the DOI | commit `4217f14` |

The commit history, the GitHub release and the Zenodo DOI are independent,
third-party-timestamped records of the date on which the benchmark folds, the
model and the analyses in this repository became public.

## What is claimed here

- the leakage-controlled benchmark: `data/processed/splits_cv.csv` and
  `data/processed/split_species.csv`, with the split-protocol audit in
  `results/leakage_ladder.csv`;
- the PalmDL model (ESM-2 view + motif/physicochemical forest) and its
  cross-validated performance;
- the cross-species transfer matrix (`results/species_matrix_auc.csv`);
- the proteome-wide FDR-controlled shortlist (`results/scan/`);
- the literature-mined audit and the structure-aware variant.

## If you build on this work

Please cite the paper and the archived release:

> PalmDL: a leakage-aware deep learning framework and benchmark for protein
> S-palmitoylation site prediction. Zenodo. https://doi.org/10.5281/zenodo.23455873

If you evaluate a new method, please use the canonical folds in
`data/processed/splits_cv.csv` (cluster-held-out) and report the protocol, so
that numbers remain comparable. The benchmark is released under CC BY 4.0; the
code under MIT (`LICENSE`, `LICENSE-DATA`).
