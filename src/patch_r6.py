"""Sixth-round cleanup: causal wording, table notes, peripheral docs, bundle README."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# ---------------- manuscript ----------------
md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")

old_dup = ("Fit and target rows exclude human background negatives; the human early-stopping split "
           "additionally contains 6,298 background negatives, unlike the other rows. Twenty percent of "
           "source clusters are held out for the diagonal and another 15% for early stopping. For the "
           "human source the early-stopping split also contained 6,298 background negatives, unlike the "
           "other rows; the 18-cell mammalian off-diagonal mean is 0.901.")
new_dup = ("Fit and target rows exclude human background negatives; the human early-stopping split "
           "additionally contains 6,298 background negatives, unlike the other rows. Twenty percent of "
           "source clusters are held out for the diagonal and another 15% for early stopping; the "
           "18-cell mammalian off-diagonal mean is 0.901.")
if old_dup in t:
    t = t.replace(old_dup, new_dup)

t = t.replace(
    "the extracellular estimate draws on three of five folds with both classes and the lumenal class (14 positives) is not interpretable",
    "the extracellular stratum contains 27 positives distributed across three folds while its pooled metrics include rows from all five folds, so it should be interpreted cautiously; the lumenal class (14 positives) is not interpretable")

t = t.replace(
    "Two design choices follow from the error analysis. First, the gap between inferred (pooled AUC 0.941) and experimental sites (0.902) is an annotation-memory effect that any user should know about.",
    "Two design choices follow from the error analysis. First, pooled AUC is 0.941 for inferred sites and 0.902 for experimentally supported sites; this evidence-tier difference matters to users and does not by itself identify its biological or annotation-related cause.")

t = t.replace(
    "The extracellular stratum pools rows from all five folds even though only three folds contain positives of both classes, so its AUC is a pooled estimate with an unbalanced fold composition. \"Known\" refers",
    "\"Known\" refers")

table6_note = ("\n**Table 6 note.** AUC and AUPR pool all outer-test predictions. Evidence strata use their "
               "positives plus all 62,842 negatives; topology strata use rows in the same compartment. The "
               "extracellular stratum has 27 positives in three folds, while its pooled metrics include rows "
               "from all five folds.\n")
anchor = "| Topology unknown | 1,801 | 0.900 | 0.512 |\n"
if anchor in t and "Table 6 note" not in t:
    t = t.replace(anchor, anchor + table6_note, 1)

md.write_text(t, encoding="utf-8")

# ---------------- root README ----------------
r = ROOT / "README.md"
x = r.read_text(encoding="utf-8")
x = x.replace("| within species (cluster-disjoint 20% holdout) | 0.900 |",
              "| within species (cluster-disjoint 20% holdout) | 0.901 |")
x = x.replace("| *A. thaliana* → mammal | 0.765 (range 0.739–0.785) |",
              "| *A. thaliana* → mammal | 0.769 (range 0.739–0.785) |")
x = x.replace("* Pooled out-of-fold AUC: experimental 0.902 vs inferred 0.941 (annotation-memory effect).",
              "* Pooled out-of-fold AUC: experimental 0.902 vs inferred 0.941; this evidence-tier difference does not identify its cause.")
x = x.replace("""* Cytoplasmic 0.858, extracellular 0.828, topology-unknown 0.900 (pooled; the extracellular
  estimate draws on three folds with both classes).""",
              """* Cytoplasmic 0.858, extracellular 0.828, topology-unknown 0.900 (pooled). The extracellular
  stratum has 27 positives in three folds while its pooled metrics use rows from all five folds.""")
x = x.replace("93% of positives are inferred", "89.7% of positives (3,969/4,424) are inferred")
x = x.replace("still `8 × 8`", "`6 × 6`")
r.write_text(x, encoding="utf-8")

# ---------------- review package README + figure legends ----------------
mp = ROOT / "src" / "make_review_package.py"
m = mp.read_text(encoding="utf-8")
m = m.replace("- 跨物种矩阵：哺乳动物间 0.883–0.948（猪→人 0.834 例外）；植物→哺乳动物 0.739–0.781",
              "- 跨物种矩阵：哺乳动物间 18 个可观察格均值 0.901、范围 0.834–0.948（猪→人 0.834 例外）；植物→哺乳动物 0.769（0.739–0.785）")
m = m.replace("- 90% 的正样本是\"by similarity\"推断注释；实验证据子集 AUC 0.907",
              "- 约 90%（3,969/4,424）的正样本是\"by similarity\"推断注释；实验证据子集 pooled AUC 0.902")
m = m.replace("- 结构特征仅在 479 个蛋白上有模型，未并入最终发布模型",
              "- 结构特征在 476 个蛋白、1,950 个位点（463 正例）上做了配对分析，未并入发布模型")
m = m.replace("蛋白级共享窗口仍概括为 `40–44%`", "蛋白级共享窗口 41.1%（随机位点级）/45.2%（蛋白级）")
m = m.replace("**Figure 2. 划分方式决定表观性能。** (a) 四种协议下的 AUC 与平衡 MCC；(b) 测试正样本中精确 31-aa 窗口出现在训练集的比例（随机/蛋白级 40–44%，聚类/物种 0）。",
              "**Figure 2. 划分方式决定表观性能。** (a) 四种协议下的 AUC 与平衡 MCC；(b) 测试正样本中精确 31-aa 窗口出现在训练集的比例（随机位点级 41.1%、蛋白级 45.2%，聚类/物种 0）。")
m = m.replace("**(b) Estimated FDR as a function of shortlist size**", "**(b) Empirical q-value as a function of shortlist size**")
mp.write_text(m, encoding="utf-8")

# ---------------- bundle README template ----------------
gb = ROOT / "src" / "make_github_bundle.py"
g = gb.read_text(encoding="utf-8")
start = g.index('README = """')
end = g.index('"""', start + 12) + 3
new_readme = '''README = """# PalmDL

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

The corrected v2.1 code and processed data live in this repository (tag v2.1). The
earlier v1.0.0 Zenodo record predates the label correction and must not be used for the
current results; a corrected Zenodo version accompanies the manuscript.

## Citation

See `CITATION.cff`; please cite the paper and the archived release.

## Licence

Code MIT (`LICENSE`), derived data CC BY 4.0 (`LICENSE-DATA`).
"""'''
g = g[:start] + new_readme + g[end:]
gb.write_text(g, encoding="utf-8")

# ---------------- final_model docstring reference ----------------
fm = ROOT / "src" / "final_model.py"
f = fm.read_text(encoding="utf-8")
f = f.replace(" (results/headline_stats.json)", "")
fm.write_text(f, encoding="utf-8")

# ---------------- reproduce note ----------------
wp = ROOT / "src" / "write_reproduce.py"
w = wp.read_text(encoding="utf-8")
w = w.replace("- `src/esm_embed.py` and `src/af_features.py` are the only steps that download data\n  (ESM-2 weights, AlphaFold models); everything else runs offline from `data/raw/`.",
              "- The initial UniProt acquisition, the literature full-text step and the model/AlphaFold\n  downloads are the only network steps; everything else runs offline from `data/raw/`.")
wp.write_text(w, encoding="utf-8")
print("sixth-round patch applied")
