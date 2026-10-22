# Reproducing the paper (v2.1)

All commands assume the repository root (`palm-dl/`) and a Python environment with
`requirements.txt` installed. On a 16-core CPU machine the pipeline takes about
three hours, dominated by ESM embedding extraction.

```bash
export PALMDL_THREADS=4          # PyTorch threads; 12 can stall on tiny batches

# Optional live refresh; overwrites the bundled 2026_03 snapshot (~1 min).
# Skip this line to reproduce the paper from the bundled tables.
# python src/download_data.py
python src/build_dataset.py      # S-palmitoyl sites + negatives (~1 min)
python src/splits.py             # cluster folds + species split (~2 min)

# main cross-validation (5 folds)
python src/run_cv.py --models pwm,lr,rf --tag cv_base
python src/run_cv.py --models cnn --tag cv_cnn
python src/esm_embed.py          # cached ESM-2 embeddings     (~20 min)
python src/run_cv.py --models esm --tag cv_esm
python src/run_cv.py --models pwm,lr,rf --subset human --tag cvh_base
python src/run_cv.py --models cnn --subset human --tag cvh_cnn
python src/run_cv.py --models esm --subset human --tag cvh_esm

# nested ensemble selection (one combination per outer fold) and splits audit
python src/ensemble_nested.py --prefix cv --out cv_ensemble
python src/ensemble_nested.py --prefix cvh --out cvh_ensemble
python src/leakage_ladder.py
python src/run_species.py --models pwm,lr,rf,cnn --tag species
python src/merge_preds.py

# tables and statistics (analysis consumes leakage_ladder and metrics_species)
python src/analysis.py
python src/stats_headline.py
python src/stats_vs_published.py

# deployed model and scan calibrator (needed by the scan and by Figure 8)
python src/final_model.py
python src/fit_rf_calibrator.py

# error structure and operating points (Figure 8 needs mutagenesis output)
python src/error_analysis.py
python src/tier_table_nested.py

# cross-species matrix v2 and the figures that consume it
python src/species_matrix_v2.py
cp results/species_matrix_v2_auc.csv results/species_matrix_auc.csv
python src/figures_extra.py      # after error_analysis: Figure 8 uses its outputs
python src/figures.py

# literature mining first: scan_analyze reads literature_sites.csv (network; ~1,000 texts)
python src/lit_mining.py --limit 1200
python src/lit_extract.py
python src/lit_eval.py

# proteome-wide scan (about one hour of CPU)
python src/scan_score.py --view rf   # forest view used in the paper
python src/scan_decoys.py
python src/scan_analyze.py
python src/fig_scan.py

# literature analyses that consume the scan output
python src/lit_recognition.py
python src/curate_literature.py  # reads all_scores and shortlists

# structure analysis (annotated human sites only; the hybrid prototype is legacy)
python src/af_features.py        # downloads AlphaFold models (~10 min with network)
python src/structure_eval.py

# documents
python src/make_supplement.py
python src/make_docx.py
```

## ESM-2 checkpoint

The language-model view needs the official fair-esm weights. Download once with
`python -c "import esm; esm.pretrained.esm2_t12_35M_UR50D()"`, then copy **both** files from
`torch.hub.get_dir()/checkpoints` into a `.cache` directory next to the `palm-dl/` checkout:

- `<workspace>/.cache/esm2_t12_35M_UR50D.pt`
- `<workspace>/.cache/esm2_t12_35M_UR50D-contact-regression.pt`

The fair-esm local loader reads the companion contact-regression file at load time even though
the downstream models do not predict contacts, so the two files must stay in the same directory.
`src/esm_embed.py` and `src/scan_embed.py` both expect this path; `src/af_features.py` downloads
AlphaFold models on demand.

## Release identity (v2.1)

Archived at Zenodo under the concept DOI **10.5281/zenodo.23455872** (version 2.1.0: 10.5281/zenodo.23537478).
The repository tag is **v2.1**; `git rev-parse v2.1` returns the exact commit used for
this manuscript. The archive uploaded to Zenodo must be built from that tag. Hashes of
the inputs and outputs quoted in the manuscript:

- `data/processed/sites.csv`  sha256 `e6a721a4e84659a246135a2f196007772178d38c19eaed1bb7744aac71d544a6`
- `data/processed/splits_cv.csv`  sha256 `19b7bd1d17052bf9c4a4e9f0f849655a78fd9948bef47d836d28e93dd87fd1ce`
- `data/processed/split_species.csv`  sha256 `c4621bfee5fdf80670c8adfdf69d3344a23d1d20e4dcdeea1023c9337d7d9963`
- `results/preds_cv_all.csv`  sha256 `b28e09c6704efe3990241b00e33ba52ff45d11a86aa2c19d8b6e0e20f2c65459`
- `results/species_matrix_auc.csv`  sha256 `f6d42da876c600fe66fcd8a2d192a7c3c5cda2b252cf106a0e1bf64a5e1735b1`
- `results/error_analysis/table6_nested.csv`  sha256 `937ba7f911c38cab414fc8d27cbcaa2f23bd5f6cd860eb86e084138af21936ba`
- `results/error_analysis/ece_nested.json`  sha256 `79b93bcf7139b1ce06d7d0f2eefbf9e4fe447f305a90eceba9b14534c8a57e93`
- `manuscript/PalmDL_manuscript.docx`  sha256 `cdadac4f350eab37b091e8d7d451b3e8dc4e882da3da05fb5037f2e0ce2326d3`
- `manuscript/manuscript.md`  sha256 `60e5f8c170a9a6aac565983c20f5ec1f3563fd0b2646eaa26d0dc3b942513235`

## Notes

- The initial UniProt acquisition, the literature full-text step and the model/AlphaFold
  downloads are the only network steps; everything else runs offline from `data/raw/`.
- Model binaries are regenerated by `src/final_model.py`; the forest-only scan
  calibrator by `src/fit_rf_calibrator.py`.
- The legacy structure scripts (`structure_hybrid.py`, `structure_rerank.py`,
  `structure_scorer.py`) belong to the superseded label set and are kept only for the
  audit trail; they are not part of the v2.1 results.
- The canonical benchmark folds are `data/processed/splits_cv.csv`; please report
  results on those folds so numbers remain comparable with the paper.
