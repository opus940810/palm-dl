"""Refresh README top sections for the S-only rebuild."""
from __future__ import annotations

from pathlib import Path

README = Path(__file__).resolve().parents[1] / "README.md"
lines = README.read_text(encoding="utf-8").split("\n")

HEAD = """# PalmDL — leakage-aware deep learning for protein S-palmitoylation site prediction

PalmDL is a two-view ensemble (ESM-2 protein-language-model embeddings + random forest on
motif/physicochemical features) for predicting whether a cysteine is S-palmitoylated, together
with a leakage-controlled benchmark built from a single, versioned UniProt release.

**Label correction (v2).** The first release filtered lipidation notes with the substring
`"palmitoyl"`, which admitted 2,381 N-/O-palmitoyl sites (35% of the positives). The corrected
dataset keeps only S-palmitoyl cysteine, contains no non-cysteine positives, and every number in
this file refers to the corrected data. The audit trail and the corrected site table are released.

## Key results

| Evaluation | n positives | AUC | MCC | Balanced MCC |
|---|---|---|---|---|
| Cluster-held-out 5-fold CV, all species | 4,424 | **0.944** (0.940–0.948) | **0.611** | 0.766 |
| Human test folds (multi-species training) | 557 | **0.914** (0.900–0.927) | 0.314 | 0.689 |
| Human-only training, human test (ESM-2 view) | 557 | 0.892 | 0.286 | 0.609 |
| Strict species holdout (random forest) | 557 | 0.885 | 0.578 | 0.614 |
| Hard negatives only (positives + same-protein Cys) | 4,424 | 0.956 | 0.678 | 0.775 |

Reference points: PWM baseline (CSS-Palm-style) 0.799 AUC; multi-kernel CNN 0.854; random forest
0.931; ESM-2 view 0.931. Published human AUC to beat: **0.900** (GPS-Palm, Brief Bioinform 2021);
PalmDL reaches 0.914, a marginal advantage whose 95% CI touches 0.900 (P(AUC ≤ 0.900) = 0.025).
Ensemble gains over either single view are +0.013 AUC (paired bootstrap, CI excludes zero).

## Theory: the code is shared across mammals and partly with plants

Single-species models evaluated on all species (`results/species_matrix_auc.csv`, Figure 7):

| transfer regime | mean AUC |
|---|---|
| within species (cluster-disjoint 20% holdout) | 0.925 |
| mammal → mammal | 0.898 |
| mammal → *A. thaliana* | 0.974 |
| *A. thaliana* → mammal | 0.774 |

Orthologue-backed transfer is nearly free (0.90), but transfer to human families with all
human-homologous clusters removed drops to 0.885. Bacterial lipidation is a different modification
(N-palmitoyl cysteine) and is excluded by construction rather than compared.

## Error structure and FDR control

* Experimental positives AUC 0.907 vs inferred 0.947 (annotation-memory effect).
* Cytoplasmic 0.887, extracellular 0.855, topology-unknown 0.898.
* FDR-controlled operating points, thresholds chosen inside each fold: 1% → precision 0.953 /
  recall 20.2%; 5% → 0.903 / 34.5%; 10% → 0.855 / 40.4%
  (`results/error_analysis/fdr_operating_points.csv`, Figure 8b).

## Interpretability and the label bug

In-silico mutagenesis (Figure 8a): cysteine at ±1 raises the score by 0.083–0.089 (cysteine
clusters), acidic residues at −1 lower it by 0.06–0.09. Mutating the central cysteine is neutral.
The same audit exposed our own filter bug: the 66 non-cysteine positives were legitimate
N6-/N-/O-palmitoyl annotations, not UniProt errors; the earlier claim is withdrawn in the manuscript.

## Split-protocol audit

The same random forest under four protocols (`results/leakage_ladder.csv`, Figure 2):

| protocol | shared test-positive windows | AUC | MCC |
|---|---|---|---|
| random site-level | 41.1% | 0.980 | 0.843 |
| protein-level | 45.2% | 0.985 | 0.810 |
| cluster-held-out | **0%** | 0.944 | 0.720 |
| species holdout | **0%** | 0.885 | 0.578 |

Cluster and species splits share no exact window at all — positives and background negatives alike.
"""

start = next(i for i, ln in enumerate(lines) if ln.startswith("## Repository layout"))
README.write_text("\n".join(HEAD.split("\n") + lines[start:]), encoding="utf-8")
print("README head replaced; total lines", len(HEAD.split("\n")) + len(lines) - start)
