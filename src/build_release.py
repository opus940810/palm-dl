"""Assemble the Zenodo-ready release bundle.

Excludes third-party full texts (copyright), AlphaFold PDB files and cached
embeddings (large, re-downloadable), and model binaries (regenerable in minutes).
Includes the exact UniProt snapshot tables so the numbers can be reproduced
offline, plus every summary table, figure, the manuscript and the service.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / "release"
VERSION = "v2.1.1"
STAGE = REL / f"palmdl-{VERSION}"

INCLUDE = [
    "src", "server", "figures", "manuscript", "deploy", "docs",
    "data/processed", "data/raw", "results/literature/review_sheet.tsv",
    "results/tables", "results/error_analysis", "results/literature/curated_sites.csv",
    "results/literature/filter_table.csv", "results/literature/metrics.json",
    "results/literature/recognition_metrics.json", "results/literature/literature_sites.csv",
    "results/scan/fdr_summary.csv", "results/scan/compartment_enrichment.csv",
    "results/scan/scan_summary.json", "results/scan/shortlist_fdr1.csv",
    "results/scan/shortlist_fdr5.csv", "results/scan/shortlist_fdr20.csv",
    "results/scan/all_scores.csv",
    "results/structure/hybrid_summary.json", "results/structure/hybrid_cv.csv",
    "results/structure/thirdview_summary.json", "results/structure/thirdview_cv.csv",
    "results/structure/paired_cv.csv", "results/structure/structure_summary.json",
    "results/structure/literature_struct_eval.json", "results/structure/scan_struct_eval.json",
    "results/structure/hybrid_apply.json", "results/structure/reranker_cv_summary.json",
    "results/preds_cv_all.csv", "results/headline_stats.json", "results/stats_vs_published.json",
    "results/leakage_ladder.csv", "results/sensitivity_residue_filter.csv",
    "results/species_matrix_auc.csv", "results/species_matrix.csv",
    "results/species_matrix_summary.json", "results/species_matrix_diagonal.csv",
    "results/species_matrix_search.json" if False else "results/ensemble_search_cvh.csv",
    "results/ensemble_search_cv.csv",
]


def copy_path(rel: str) -> int:
    src = ROOT / rel
    dst = STAGE / rel
    if not src.exists():
        print("  missing (skipped):", rel)
        return 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
        return sum(1 for _ in dst.rglob("*") if _.is_file())
    shutil.copy2(src, dst)
    return 1


def main() -> None:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    total = 0
    for rel in INCLUDE:
        total += copy_path(rel)
    print("files staged:", total)

    (STAGE / "LICENSE").write_text(
        "MIT License\n\nCopyright (c) 2026 The PalmDL authors\n\nPermission is hereby granted, free "
        "of charge, to any person obtaining a copy of this software and associated documentation "
        "files (the \"Software\"), to deal in the Software without restriction, including without "
        "limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or "
        "sell copies of the Software, and to permit persons to whom the Software is furnished to do "
        "so, subject to the following conditions:\n\nThe above copyright notice and this permission "
        "notice shall be included in all copies or substantial portions of the Software.\n\nTHE "
        "SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.\n",
        encoding="utf-8")
    (STAGE / "LICENSE-DATA").write_text(
        "Derived data in this archive (site tables, fold assignments, predictions, scores, figures) "
        "are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0).\n\n"
        "Source data: UniProt release 2026_03, CC BY 4.0 (https://www.uniprot.org). AlphaFold "
        "structures are not redistributed here; they are available from the AlphaFold DB under "
        "CC BY 4.0.\n",
        encoding="utf-8")
    (STAGE / "requirements.txt").write_text(
        "numpy>=2.0\npandas>=2.2\ntorch==2.9.0\nscikit-learn>=1.5\nscipy>=1.13\n"
        "matplotlib>=3.8\nbiopython>=1.83\nfair-esm==2.0.0\npython-docx>=1.1\n",
        encoding="utf-8")
    (STAGE / "CITATION.cff").write_text(
        "cff-version: 1.2.0\nmessage: If you use PalmDL, please cite the paper and this archive.\n"
        "title: \"PalmDL: a leakage-aware deep learning framework and benchmark for protein "
        "S-palmitoylation site prediction\"\nversion: " + VERSION + "\n"
        "authors:\n  - family-names: \"[to be completed]\"\n    given-names: \"[to be completed]\"\n"
        "license: MIT\nrepository-code: \"https://doi.org/[ZENODO-DOI]\"\n",
        encoding="utf-8")
    (STAGE / ".zenodo.json").write_text(json.dumps({
        "title": "PalmDL: leakage-aware deep learning and benchmark for protein S-palmitoylation site prediction",
        "upload_type": "software",
        "description": ("Code, processed data, fold assignments and results for PalmDL. Includes a "
                        "leakage-controlled benchmark (cluster- and species-held-out folds), the "
                        "ESM-2 + random-forest ensemble, a structure-aware variant, a "
                        "proteome-wide FDR-controlled shortlist, and a prediction service."),
        "creators": [{"name": "[Family name], [Given name]", "affiliation": "[Affiliation]"}],
        "license": "MIT",
        "keywords": ["S-palmitoylation", "post-translational modification", "deep learning",
                     "protein language model", "benchmark leakage", "cross-species"],
        "related_identifiers": [{"relation": "isSupplementTo",
                                 "identifier": "10.xxxx/[manuscript DOI]",
                                 "resource_type": "publication-article"}],
        "communities": [], "access_right": "open",
    }, indent=2), encoding="utf-8")

    (STAGE / "README.md").write_text(
        "# PalmDL release " + VERSION + "\n\n"
        "Code, processed data, fold assignments, results and the prediction service for\n"
        "*PalmDL: a leakage-aware deep learning framework and benchmark for protein\n"
        "S-palmitoylation site prediction*.\n\n"
        "## Layout\n\n"
        "- `src/` - 25 numbered pipeline scripts (download, dataset, splits, models, scan, literature, structure)\n"
        "- `server/` - standard-library prediction service (HTTP + batch; `--structure` switch)\n"
        "- `data/processed/` - labelled sites and **the canonical fold assignments**\n"
        "- `data/raw/` - the exact UniProt 2026_03 snapshot used (download metadata included)\n"
        "- `results/` - per-fold predictions, summary tables, proteome-wide scores, shortlists\n"
        "- `figures/`, `manuscript/` - main text and supplementary material\n"
        "- `deploy/` - Docker/compose/nginx recipes for hosting the service\n\n"
        "## Quick start\n\n"
        "```bash\npip install -r requirements.txt\n"
        "python server/app.py --batch server/example.fa --out example.csv --threshold 0.5\n"
        "python server/app.py                # http://127.0.0.1:8765\n```\n\n"
        "## Reproducing the paper\n\n"
        "Run the scripts in numeric order (see `docs/REPRODUCE.md`). Model binaries are not bundled "
        "(size); they are regenerated by `src/final_model.py` (~3 min) and `src/structure_hybrid.py` "
        "(~10 min) after downloading the ESM-2 35M weights with fair-esm.\n\n"
        "## Licence\n\nMIT for code; CC BY 4.0 for derived data. UniProt 2026_03 is CC BY 4.0.\n",
        encoding="utf-8")
    (STAGE / "docs").mkdir(exist_ok=True)

    # checksums
    lines = []
    for p in sorted(STAGE.rglob("*")):
        if p.is_file() and p.name != "checksums.sha256":
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            lines.append(f"{h}  {p.relative_to(STAGE).as_posix()}")
    (STAGE / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    archive = REL / f"palmdl-{VERSION}.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(STAGE.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(STAGE).as_posix())
    size = archive.stat().st_size / 1e6
    print(f"archive: {archive} ({size:.1f} MB, {len(lines)} files)")


if __name__ == "__main__":
    main()
