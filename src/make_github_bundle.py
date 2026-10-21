"""Build a lean, push-ready GitHub repository tree (dist/github).

Big generated artifacts (out-of-fold predictions, proteome-wide score table,
UniProt raw snapshot) are deliberately excluded: they live in the Zenodo archive.
The repository carries code, docs, figures, manuscript, the canonical folds and
all summary tables, so it stays well under GitHub's file-size limits.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "dist" / "github"

INCLUDE = [
    "src", "server", "deploy", "docs", "figures", "manuscript",
    "data/processed",
    "results/tables", "results/error_analysis",
    "results/literature/curated_sites.csv", "results/literature/filter_table.csv",
    "results/literature/metrics.json", "results/literature/recognition_metrics.json",
    "results/literature/literature_sites.csv", "results/literature/curated_overlap.json",
    "results/scan/fdr_summary.csv", "results/scan/compartment_enrichment.csv",
    "results/scan/scan_summary.json", "results/scan/shortlist_fdr1.csv",
    "results/scan/shortlist_fdr5.csv", "results/scan/shortlist_fdr20.csv",
    "results/structure/hybrid_summary.json", "results/structure/hybrid_cv.csv",
    "results/structure/thirdview_summary.json", "results/structure/thirdview_cv.csv",
    "results/structure/paired_cv.csv", "results/structure/structure_summary.json",
    "results/structure/literature_struct_eval.json", "results/structure/scan_struct_eval.json",
    "results/structure/hybrid_apply.json", "results/structure/reranker_cv_summary.json",
    "results/headline_stats.json", "results/stats_vs_published.json",
    "results/leakage_ladder.csv", "results/sensitivity_residue_filter.csv",
    "results/species_matrix_auc.csv", "results/species_matrix.csv",
    "results/species_matrix_summary.json", "results/species_matrix_diagonal.csv",
    "results/ensemble_search_cv.csv", "results/ensemble_search_cvh.csv",
]

GITIGNORE = """# Python
__pycache__/
*.py[cod]
.venv/
venv/

# model binaries and caches (regenerate with src/final_model.py, src/structure_hybrid.py)
models/
.cache/
.pylibs/

# large generated artifacts (available in the Zenodo archive)
results/preds_cv_all.csv
results/scan/all_scores.csv
results/scan/chunk_*.csv
results/literature/xml/
results/structure/af/
data/raw/*.tsv

# office / OS noise
~$*
.DS_Store
Thumbs.db
"""

GITATTRIBUTES = """* text=auto
*.png binary
*.docx binary
*.joblib binary
*.pt binary
"""

README = """# PalmDL

[![DOI](https://zenodo.org/badge/DOI/[ZENODO-DOI].svg)](https://doi.org/[ZENODO-DOI])
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

The corrected v2.1 code and processed data live in this repository (tag v2.1) and are
archived at Zenodo under the concept DOI 10.5281/zenodo.23455872 (version 2.1.0: 10.5281/zenodo.23537476). The earlier
v1.0.0 record (10.5281/zenodo.23455873) predates the label correction and must not be used
for the current results.

## Citation

See `CITATION.cff`; please cite the paper and the archived release.

## Licence

Code MIT (`LICENSE`), derived data CC BY 4.0 (`LICENSE-DATA`).
"""

MIT = """MIT License

Copyright (c) 2026 The PalmDL authors

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
"""

CCBY = """Derived data in this repository (site tables, fold assignments, summary
predictions, scores, figures) are released under the Creative Commons Attribution
4.0 International licence (CC BY 4.0).

Source data: UniProt release 2026_03 (CC BY 4.0, https://www.uniprot.org).
AlphaFold structures are not redistributed here; they are available from the
AlphaFold DB under CC BY 4.0.
"""

CITATION = """cff-version: 1.2.0
message: If you use PalmDL, please cite the paper and the archived release.
title: "PalmDL: a leakage-aware deep learning framework and benchmark for protein S-palmitoylation site prediction"
version: 1.0.0
authors:
  - family-names: "[to be completed]"
    given-names: "[to be completed]"
license: MIT
repository-code: "https://github.com/[ORG]/palm-dl"
identifiers:
  - type: doi
    value: "[ZENODO-DOI]"
"""

REQUIREMENTS = """numpy>=2.0
pandas>=2.2
torch==2.9.0
scikit-learn>=1.5
scipy>=1.13
matplotlib>=3.8
biopython>=1.83
fair-esm==2.0.0
python-docx>=1.1
"""

ZENODO = {
    "title": "PalmDL: leakage-aware deep learning and benchmark for protein S-palmitoylation site prediction",
    "upload_type": "software",
    "description": ("Code, processed data, fold assignments and results for PalmDL: a "
                    "leakage-controlled benchmark, an ESM-2 + random-forest ensemble, a "
                    "structure-aware variant, a proteome-wide FDR-controlled shortlist and a "
                    "prediction service."),
    "creators": [{"name": "[Family name], [Given name]", "affiliation": "[Affiliation]"}],
    "license": "MIT",
    "keywords": ["S-palmitoylation", "post-translational modification", "deep learning",
                 "benchmark leakage"],
    "access_right": "open",
}



CHANGELOG = """# Changelog

## v2 (19 October 2026) - label correction

- The v1 dataset filtered lipidation notes with the substring `"palmitoyl"`, admitting
  2,381 N-/O-palmitoyl sites (35% of positives). v2 keeps only S-palmitoyl cysteine on a
  cysteine residue: **4,424 positives** (449 experimental, 3,975 inferred), 67,266 sites,
  21,568 proteins, 786 organisms, 557 human positives.
- Splits: every protein (including background negatives) is clustered; groups sharing an
  exact window are merged. Cross-fold identical windows: **0** (v1: 751 windows).
- Results corrected throughout: cluster CV AUC 0.944 (0.940-0.948), MCC 0.611; human
  0.914 (0.900-0.927); species holdout 0.885; hard negatives 0.956.
- The "UniProt annotation error" claim for the 66 non-cysteine sites is withdrawn: those
  were legitimate N6-/N-/O-palmitoyl annotations admitted by our filter.
- Cross-species matrix reduced to six species; bacterial N-palmitoyl sites removed.
- Proteome scan rescored with the forest view (ESM view not tractable at proteome scale
  on the available CPU host); tiers reported as candidate rankings.
- ESM-2 view omitted from the scan; the released prediction service still uses the full
  sequence ensemble.
"""


def main() -> None:
    if OUT.exists():
        # keep .git so the published history survives regeneration
        for child in OUT.iterdir():
            if child.name == ".git":
                continue
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    OUT.mkdir(parents=True, exist_ok=True)
    n = 0
    for rel in INCLUDE:
        src, dst = ROOT / rel, OUT / rel
        if not src.exists():
            print("  missing (skipped):", rel)
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dst)
        n += 1
    (OUT / "LICENSE").write_text(MIT, encoding="utf-8")
    (OUT / "LICENSE-DATA").write_text(CCBY, encoding="utf-8")
    (OUT / "CITATION.cff").write_text(CITATION, encoding="utf-8")
    (OUT / "requirements.txt").write_text(REQUIREMENTS, encoding="utf-8")
    (OUT / ".zenodo.json").write_text(json.dumps(ZENODO, indent=2), encoding="utf-8")
    (OUT / "CHANGELOG.md").write_text(CHANGELOG, encoding="utf-8")
    (OUT / ".gitignore").write_text(GITIGNORE, encoding="utf-8")
    (OUT / ".gitattributes").write_text(GITATTRIBUTES, encoding="utf-8")
    (OUT / "README.md").write_text(README, encoding="utf-8")
    size = sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file()) / 1e6
    files = sum(1 for p in OUT.rglob("*") if p.is_file())
    biggest = sorted(((p.stat().st_size, p) for p in OUT.rglob("*") if p.is_file()), reverse=True)[:5]
    print(f"github tree: {files} files, {size:.1f} MB")
    for s, p in biggest:
        print(f"  {p.relative_to(OUT)}  {s / 1e6:.1f} MB")
    print("next: push dist/github to a new repository (see docs/GITHUB_RELEASE.md)")


if __name__ == "__main__":
    main()
