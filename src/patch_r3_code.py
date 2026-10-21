"""Final third-round code/doc fixes: structure paragraph, 4.8 heading, server UI, reproduce."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# manuscript
md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")
old = t[t.index("It does not cure the domain shift of §2.10."):]
old = old[:old.index("\n\n")]
t = t.replace(old, "No structural model was run for the literature sites in this revision, so whether structure would help there is untested. The earlier claim that structure cannot close the domain shift is therefore withdrawn: the paired comparison only establishes a modest gain where annotated human sites have AlphaFold models.")
t = t.replace("### 4.8 FDR-controlled operating points", "### 4.8 Empirical validation-target operating points")
md.write_text(t, encoding="utf-8")

# server app: remove structure UI and stale claims
app = ROOT / "server" / "app.py"
s = app.read_text(encoding="utf-8")
s = s.replace('''<label>threshold <input type="number" name="threshold" step="0.01" value="0.5"></label>
&nbsp;<label><input type="checkbox" name="structure" value="1"> use AlphaFold structure if available</label>
&nbsp;<button type="submit">Predict</button>''', '''<label>threshold <input type="number" name="threshold" step="0.01" value="0.5"></label>
&nbsp;<button type="submit">Predict</button>''')
s = s.replace('''<p class="hint">Paste one or more protein sequences in FASTA format. Every cysteine is scored;
the table is ranked by probability. Model: ESM-2 embedding view + motif/physicochemical
random forest, Platt-calibrated. Default threshold is the natural-prevalence operating point;
use 0.957 for the proteome-wide 5% FDR tier. If the record name is a UniProt accession and a
structure is available, the structure-aware path is used for modelled cysteines.</p>''', '''<p class="hint">Paste one or more protein sequences in FASTA format. Every cysteine is scored and
the table is ranked by probability. The service uses the fixed RF+ESM logit-average model with a
Platt calibration; it is not the nested cross-validation benchmark. The natural-prevalence
operating point is 0.36; higher thresholds are more conservative. Proteome-scan tiers in the
manuscript are on a different (forest-only) scale and are not comparable with these probabilities.</p>''')
s = s.replace("Modes\n  sequence (default)  RF + ESM ensemble, available for any FASTA input.\n  structure           if an AlphaFold model exists for the record name\n                      (a UniProt accession), the hybrid structure-aware RF is\n                      used for modelled cysteines and the sequence path is used\n                      for the rest.  Models are taken from a local cache and,\n                      with --fetch, downloaded on demand.",
              "Mode\n  sequence  RF + ESM logit-average ensemble with Platt calibration, available for any\n            FASTA input. The legacy structure path is retired: requests cannot enable it.")
s = s.replace("structure=models rf_hybrid_struct", "structure=models rf_hybrid_struct")  # no-op
app.write_text(s, encoding="utf-8")

# server README: drop residual FDR/structure hints
sr = ROOT / "server" / "README.md"
r = sr.read_text(encoding="utf-8")
r = r.replace("""- `0.36` — natural-prevalence operating point of the released model.
- `0.5` — conservative default for ranking.
- The proteome-scan tiers in the manuscript (Table 8) are on the **forest-view** scale and are
  not comparable with service probabilities; use them only to rank candidates inside that table.""",
"""- `0.36` — natural-prevalence operating point of the deployed logit-average model.
- `0.5` — conservative default for ranking.
- The nested cross-validation benchmark and the deployed model are different algorithms: the
  benchmark selects members per fold and averages within-fold ranks, whereas the service fixes
  RF+ESM and averages logits. Proteome-scan tiers are on a forest-only scale and are not
  comparable with service probabilities.""")
sr.write_text(r, encoding="utf-8")

# reproduce doc
rp = ROOT / "docs" / "REPRODUCE.md"
d = rp.read_text(encoding="utf-8")
d = d.replace("""# ensembles, tables, statistics, figures
python src/ensemble.py --members rf,esm --out cv_ensemble
python src/ensemble.py --members rf,esm --out cvh_ensemble""",
"""# nested ensemble selection (one combination per outer fold)
python src/run_cv.py --models pwm,lr,rf --subset human --tag cvh_base
python src/ensemble_nested.py --prefix cv --out cv_ensemble
python src/ensemble_nested.py --prefix cvh --out cvh_ensemble""")
d = d.replace("""# cross-species matrix
python src/species_matrix.py
python src/diag_fix.py           # honest diagonal (cluster holdout)
python src/figures_extra.py""",
"""# cross-species matrix (v2: frozen targets, fit/val cluster exclusion)
python src/species_matrix_v2.py
cp results/species_matrix_v2_auc.csv results/species_matrix_auc.csv
python src/figures_extra.py""")
d = d.replace("""python src/error_analysis.py""", """python src/final_model.py
python src/fit_rf_calibrator.py
python src/error_analysis.py      # Figure 8 depends on this
python src/tier_table_nested.py""")
rp.write_text(d, encoding="utf-8")
print("third-round code/doc fixes applied")
