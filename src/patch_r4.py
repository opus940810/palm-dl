"""Fourth-round cleanup: README entry, ordering, text residuals, server threshold."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# ---------- manuscript ----------
md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")


def rep(old: str, new: str) -> None:
    global t
    if old not in t:
        print("MISS:", old[:70])
        return
    t = t.replace(old, new)


rep("The plant model reaches 0.74–0.78 on mammals", "The plant model reaches 0.739–0.785 on mammals")
rep("mammal→plant 0.964–0.984", "mammal→plant 0.964–0.984")
rep("plant model reaches 0.74–0.78 on mammals", "plant model reaches 0.739–0.785 on mammals")
rep("### 2.12 A structure-aware prototype: useful inside the annotation distribution, not outside it",
    "### 2.12 A paired structural analysis on annotated human sites")
rep("Two design choices follow from the error analysis. First, the gap between inferred (AUC 0.947) and experimental sites (0.907)",
    "Two design choices follow from the error analysis. First, the gap between inferred (pooled AUC 0.941) and experimental sites (0.902)")
rep("With orthologues available, most mammalian pairs transfer at 0.83–0.95 (18-cell mean 0.901)",
    "Across the sampled species, most mammalian pairs transfer at 0.83–0.95 (18-cell mean 0.901)")
rep("One model per source species is trained on its positives and same-protein negatives (no human background); 20% of source clusters are held out for the diagonal and another 15% for early stopping.",
    "Fit and target rows exclude human background negatives; the human early-stopping split additionally contains 6,298 background negatives, unlike the other rows. Twenty percent of source clusters are held out for the diagonal and another 15% for early stopping.")
rep("### 4.7 Cross-species transfer matrix", "### 4.7 Cross-species transfer matrix")
rep("that species' positives and its same-protein negatives; human background negatives were excluded.",
    "that species' positives and its same-protein negatives; fit and target rows exclude human background negatives, while the human early-stopping split also contains 6,298 background negatives (disclosed in §2.7).")
rep("FDR is estimated by target–decoy comparison against 20,000 shuffled-flank decoys scored with the forest view; the decoy tail is thin, so these tiers are candidate rankings rather than proven error rates.",
    "Empirical q-values are estimated by target–decoy comparison against 20,000 shuffled-flank decoys scored with the forest view; the decoy tail is thin, so these tiers are candidate rankings rather than proven error rates. Known-site recovery uses the 557 annotated human sites as its denominator.")
rep("(b) Estimated FDR as a function of shortlist size, with the five reported tiers marked.",
    "(b) Empirical q-value as a function of shortlist size, with the five reported tiers marked.")

# remove duplicated old matrix paragraph (line 106)
lines = t.split("\n")
out = []
for ln in lines:
    if ln.startswith("Diagonal entries are from a protein-cluster-disjoint 20% holdout of the same species; off-diagonal entries use a model trained on all source rows"):
        continue
    out.append(ln)
t = "\n".join(out)
md.write_text(t, encoding="utf-8")

# ---------- root README ----------
r = ROOT / "README.md"
x = r.read_text(encoding="utf-8")
old_block = """& $py palm-dl\\src\\run_species.py --models pwm,lr,rf,cnn --tag species
& $py palm-dl\\src\\ensemble.py --members rf,esm --out cv_ensemble
& $py palm-dl\\src\\ensemble.py --members rf,esm --out cvh_ensemble
& $py palm-dl\\src\\merge_preds.py
& $py palm-dl\\src\\analysis.py
& $py palm-dl\\src\\stats_headline.py
& $py palm-dl\\src\\figures.py
& $py palm-dl\\src\\make_docx.py"""
new_block = """& $py palm-dl\\src\\run_species.py --models pwm,lr,rf,cnn --tag species
& $py palm-dl\\src\\run_cv.py --models pwm,lr,rf --subset human --tag cvh_base
& $py palm-dl\\src\\ensemble_nested.py --prefix cv --out cv_ensemble
& $py palm-dl\\src\\ensemble_nested.py --prefix cvh --out cvh_ensemble
& $py palm-dl\\src\\leakage_ladder.py
& $py palm-dl\\src\\merge_preds.py
& $py palm-dl\\src\\analysis.py
& $py palm-dl\\src\\stats_headline.py
& $py palm-dl\\src\\error_analysis.py
& $py palm-dl\\src\\figures.py
& $py palm-dl\\src\\make_docx.py

The full ordered pipeline, including the matrix v2, scan and literature steps, is in
`docs/REPRODUCE.md`; use that document rather than this abbreviated block."""
if old_block in x:
    x = x.replace(old_block, new_block)
else:
    print("README block MISS")
x = x.replace("4,424 positive sites (449 experimental, 3,975 inferred)", "4,424 positive sites (455 experimental, 3,969 inferred)")
x = x.replace("""- **Ensemble selection** (all 31 member subsets) uses pooled validation AUC only; test data are
  never used for model selection. `rf+esm` was selected and is also the best under human-only CV.
- **Thresholds** come from an inner, cluster-disjoint validation split (MCC-optimal), separately
  for natural and class-balanced evaluation. Ensemble scores are Platt-calibrated (ECE 0.040).""",
"""- **Ensemble selection** is nested per outer fold: the 31 subsets are ranked on that fold's
  validation rows only, and the winner is calibrated and applied inside the same fold. Four folds
  select RF+ESM, one selects LR+RF+CNN+ESM.
- **Thresholds** come from an inner, cluster-disjoint validation split (MCC-optimal), separately
  for natural and class-balanced evaluation. Ensemble scores are Platt-calibrated (mean fold ECE
  0.024; pooled 0.018).""")
r.write_text(x, encoding="utf-8")

# ---------- server ----------
app = ROOT / "server" / "app.py"
s = app.read_text(encoding="utf-8")
s = s.replace('STRICT_THRESHOLD = 0.957  # 5% FDR at proteome scale (manuscript Table 8)\n', "")
s = s.replace("""            cls = ' class="hi"' if s["probability"] >= STRICT_THRESHOLD else ""\n""", "")
s = s.replace("""            rows.append(f"<tr{cls}><td>{rec['protein']}</td>""", """            rows.append(f"<tr><td>{rec['protein']}</td>""")
app.write_text(s, encoding="utf-8")

sr = ROOT / "server" / "README.md"
v = sr.read_text(encoding="utf-8")
v = v.replace("""- The cross-validated model used within-fold **rank** averaging, while this service averages view
  **logits** (needed for single-sequence scoring). On the same folds the logit average scores AUC
  0.940 versus 0.944 for the rank average under the global combination, so service probabilities are
  slightly less sharp than the cross-validated ensemble.""",
"""- The nested benchmark chooses members per fold and averages within-fold **ranks**; this service
  fixes RF+ESM and averages **logits**. For the fixed two-view combination the logit average scored
  AUC 0.940 versus 0.944 for the rank average in a sensitivity check; this does not establish that
  the deployed probabilities are calibrated like the benchmark.
- The structure path is retired and the interface no longer exposes it; the AlphaFold analysis is
  an annotated-human paired experiment reported in the manuscript, not a deployed feature.""")
sr.write_text(v, encoding="utf-8")

# final model docstring
fm = ROOT / "src" / "final_model.py"
f = fm.read_text(encoding="utf-8")
f = f.replace("difference < .003", "difference about .004 (fixed two-view sensitivity)")
f = f.replace("results/ensemble_logit_check.json", "results/headline_stats.json")
fm.write_text(f, encoding="utf-8")
print("fourth-round cleanup applied")
