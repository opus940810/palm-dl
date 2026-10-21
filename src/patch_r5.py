"""Fifth-round sync: README/server docs, final_model docstring, reproduce order."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# ---------------- manuscript: extracellular pooled note + Ref12 ----------------
md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")
t = t.replace(
    "\"Known\" refers to the 557 S-palmitoyl cysteine sites annotated in UniProt 2026_03 for human proteins;",
    "The extracellular stratum pools rows from all five folds even though only three folds contain positives of both classes, so its AUC is a pooled estimate with an unbalanced fold composition. \"Known\" refers to the 557 S-palmitoyl cysteine sites annotated in UniProt 2026_03 for human proteins;")
md.write_text(t, encoding="utf-8")

# ---------------- root README ----------------
r = ROOT / "README.md"
x = r.read_text(encoding="utf-8")
x = x.replace("## Theory: the code is shared across mammals and partly with plants",
              "## Theory: cross-species transfer is heterogeneous")
x = x.replace("""| within species (cluster-disjoint 20% holdout) | 0.900 |
| mammal → mammal | 0.896 (range 0.834–0.948) |
| mammal → *A. thaliana* | 0.978 |
| *A. thaliana* → mammal | 0.765 (range 0.739–0.781) |

Orthologue-backed transfer is nearly free (0.90), but transfer to human families with all
human-homologous clusters removed drops to 0.885. Bacterial lipidation is a different modification
(N-palmitoyl cysteine) and is excluded by construction rather than compared.""",
"""| within species (cluster-disjoint 20% holdout) | 0.900 |
| mammal → mammal (18 observable off-diagonal cells) | 0.901 (range 0.834–0.948) |
| mammal → *A. thaliana* | 0.978 |
| *A. thaliana* → mammal | 0.765 (range 0.739–0.785) |

Transfer is heterogeneous: Pig→Human is 0.834, and four cells have fewer than 20 eligible test
positives after cluster removal and are left empty. The design does not isolate the causal
contribution of orthologue coverage, and bacterial lipidation (N-palmitoyl cysteine) is excluded
rather than compared.""")
x = x.replace("""* Experimental positives AUC 0.907 vs inferred 0.947 (annotation-memory effect).
* Cytoplasmic 0.887, extracellular 0.855, topology-unknown 0.898.""",
"""* Pooled out-of-fold AUC: experimental 0.902 vs inferred 0.941 (annotation-memory effect).
* Cytoplasmic 0.858, extracellular 0.828, topology-unknown 0.900 (pooled; the extracellular
  estimate draws on three folds with both classes).""")
x = x.replace("## Error structure and FDR control", "## Error structure and validation-target operating points")
x = x.replace("93% of positives are inferred", "90% of positives are inferred")
x = x.replace("10. **`src/ensemble.py`**", "10. **`src/ensemble_nested.py`**")
r.write_text(x, encoding="utf-8")

# pipeline table in README
lines = x.split("\n")
for i, ln in enumerate(lines):
    if "ensemble.py" in ln and "|" in ln:
        lines[i] = ln.replace("ensemble.py", "ensemble_nested.py").replace("`", "`")
    if "species_matrix.py" in ln and "|" in ln:
        lines[i] = ln.replace("species_matrix.py", "species_matrix_v2.py")
(r).write_text("\n".join(lines), encoding="utf-8")

# ---------------- server README (full rewrite of stale parts) ----------------
sr = ROOT / "server" / "README.md"
sr.write_text("""# PalmDL prediction service

Standard-library HTTP service that scores every cysteine of one or more FASTA proteins with the
released PalmDL model (ESM-2 view + motif/physicochemical random forest, Platt-calibrated). The
service is **sequence-only**; the legacy structure path is retired.

## Run

```powershell
$py = "<python.exe>"
$env:PYTHONPATH = "<workspace>\\.pylibs"
$env:PALMDL_THREADS = "4"

# interactive server
& $py palm-dl\\server\\app.py                 # http://127.0.0.1:8765
& $py palm-dl\\server\\app.py --port 9000     # custom port

# batch mode (no server)
& $py palm-dl\\server\\app.py --batch palm-dl\\server\\example.fa `
      --out palm-dl\\server\\example_predictions.csv --threshold 0.36
```

The first request loads the model (a few seconds); subsequent requests reuse it.

## Endpoints

| Method | Path | Body | Response |
|---|---|---|---|
| GET | `/` | — | HTML form |
| GET | `/health` | — | model metadata |
| POST | `/predict` | form fields `fasta`, `threshold`; or JSON `{"fasta": "...", "threshold": 0.36}` | ranked table (HTML) or JSON |

Example:

```powershell
$body = @{ fasta = (Get-Content palm-dl\\server\\example.fa -Raw); threshold = 0.36 } | ConvertTo-Json
Invoke-RestMethod -Uri http://127.0.0.1:8765/predict -Method Post -Body $body -ContentType "application/json"
```

## Structure

The v1 service exposed an optional structure-aware path backed by a hybrid forest trained on the
superseded label set. v2 closes that path in both the CLI and HTTP handlers; `--structure` is
accepted for backward compatibility but ignored. The AlphaFold analysis in the paper is a paired
experiment on annotated human sites (forest AUC 0.890 → 0.911), not a deployed feature.

## Choosing a threshold

- `0.36` — natural-prevalence operating point of the deployed logit-average model.
- `0.5` — conservative default for ranking.
- The nested cross-validation benchmark and the deployed model are different algorithms: the
  benchmark selects members per fold and averages within-fold ranks, whereas the service fixes
  RF+ESM and averages logits. Proteome-scan tiers are on a forest-only scale and are not
  comparable with service probabilities.

## Notes

- For the fixed two-view combination the logit average scored AUC 0.940 versus 0.944 for the rank
  average in a sensitivity check (difference about 0.004); this does not establish that deployed
  probabilities are calibrated like the benchmark, whose primary estimate is the nested 0.938.
- No authentication and no upload limit beyond memory: intended for local or intranet use. Put it
  behind a reverse proxy if exposed.
""", encoding="utf-8")

# ---------------- final_model docstring ----------------
fm = ROOT / "src" / "final_model.py"
f = fm.read_text(encoding="utf-8")
f = f.replace("changes AUC by <0.003", "changes AUC by about 0.004 in a fixed two-view sensitivity check")
fm.write_text(f, encoding="utf-8")

# ---------------- review package figure legends ----------------
mp = ROOT / "src" / "make_review_package.py"
m = mp.read_text(encoding="utf-8")
m = m.replace("**Figure 1. 数据构成。** (a) 贡献位点最多的物种（人/小鼠高亮）；(b) 证据层级（477 实验 vs 6,328 推断，对数坐标）；(c) 正位点的拓扑注释。",
              "**Figure 1. 数据构成。** (a) 贡献位点最多的物种（人/小鼠高亮）；(b) 证据层级（455 实验 vs 3,969 推断，对数坐标）；(c) 正位点的拓扑注释。")
m = m.replace("**Figure 7. 棕榈酰化\"密码\"是真核共有的。** 单物种模型的 8×8 迁移矩阵",
              "**Figure 7. 跨物种迁移是异质的。** 单物种模型的 6×6 迁移矩阵")
m = m.replace("对角线来自同物种 20% 聚类留出。", "对角线来自同物种 20% 聚类留出；四个空格因可用正例不足 20 而留空。")
mp.write_text(m, encoding="utf-8")
print("fifth-round docs patched")
