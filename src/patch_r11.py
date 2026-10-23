"""Tenth-round P3 sync: template strings and archive/version labelling."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = "10.5281/zenodo.23550005"
DOCREV = "10.5281/zenodo.23570441"

# 1) review-package data dictionary 6x6
mp = ROOT / "src" / "make_review_package.py"
m = mp.read_text(encoding="utf-8")
m = m.replace("| `species_matrix_auc.csv` | 8×8 跨物种 AUC 矩阵 |",
              "| `species_matrix_auc.csv` | 6×6 跨物种 AUC 矩阵 |")
# review README archive labelling
m = m.replace("**代码与归档 / Code and archive:** https://github.com/opus940810/palm-dl （tag v2.1）｜ Zenodo concept DOI: 10.5281/zenodo.23455872（最新版 2.1.2）",
              f"**代码与归档 / Code and archive:** https://github.com/opus940810/palm-dl （tag v2.1）｜ "
              f"冻结实验归档 v2.1.2: {FROZEN}｜文档修订版 v2.1.3: {DOCREV}")
mp.write_text(m, encoding="utf-8")

# 2) bundle metadata template
gb = ROOT / "src" / "make_github_bundle.py"
g = gb.read_text(encoding="utf-8")
g = g.replace('"structure-aware variant, a proteome-wide FDR-controlled shortlist and a "',
              '"six-species transfer matrix, a proteome-wide empirical candidate list, and a "')
g = g.replace('"prediction service."', '"sequence-only prediction service."')
gb.write_text(g, encoding="utf-8")

# 3) REPRODUCE labelling
wp = ROOT / "src" / "write_reproduce.py"
w = wp.read_text(encoding="utf-8")
w = w.replace('Archived at Zenodo under the concept DOI **{CONCEPT}**; the exact snapshot is pinned by tag v2.1.',
              f"Archived at Zenodo under the concept DOI **{{CONCEPT}}**. The frozen experiment archive is "
              f"**v2.1.2** ({FROZEN}); **v2.1.3** ({DOCREV}) is a documentation revision of the same results. "
              f"The working tree is pinned by GitHub tag v2.1.")
wp.write_text(w, encoding="utf-8")
print("tenth-round P3 sync applied")
