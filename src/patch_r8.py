"""Eighth-round fixes: public DOI template, ESM path, raw in git, reproduce refresh."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# 1) bundle template: correct version DOI + include raw snapshot
gb = ROOT / "src" / "make_github_bundle.py"
g = gb.read_text(encoding="utf-8")
g = g.replace("10.5281/zenodo.23537476", "10.5281/zenodo.23537478")
g = g.replace('"src", "server", "deploy", "docs", "figures", "manuscript",\n    "data/processed",',
              '"src", "server", "deploy", "docs", "figures", "manuscript",\n    "data/processed", "data/raw",')
g = g.replace("data/raw/*.tsv\n", "")
# README wording: raw snapshot is now in the repository
g = g.replace("The exact UniProt 2026_03 snapshot used for the paper is included under `data/raw/`; re-running",
              "The exact UniProt 2026_03 snapshot used for the paper is included under `data/raw/` (and in the\n"
              "archived release); re-running")
gb.write_text(g, encoding="utf-8")

# 2) root README pipeline row
r = ROOT / "README.md"
lines = r.read_text(encoding="utf-8").split("\n")
for i, ln in enumerate(lines):
    if "species_matrix_v2.py" in ln and "diag_fix" in ln:
        lines[i] = "| 14 | `src/species_matrix_v2.py` | 6 × 6 cross-species transfer matrix (v2) |"
r.write_text("\n".join(lines), encoding="utf-8")

# 3) Chinese figure 2 caption wording
c = ROOT / "中文说明.md"
s = c.read_text(encoding="utf-8")
s = s.replace("随机/蛋白级 40–44%", "随机位点级 41.1%、蛋白级 45.2%")
s = s.replace("随机/蛋白级 41–45%", "随机位点级 41.1%、蛋白级 45.2%")
s = s.replace("（随机/蛋白级 40–44%，聚类/物种 0）", "（随机位点级 41.1%、蛋白级 45.2%，聚类/物种 0）")
c.write_text(s, encoding="utf-8")

# 4) reproduce: ESM path + optional live refresh
wp = ROOT / "src" / "write_reproduce.py"
w = wp.read_text(encoding="utf-8")
w = w.replace("python src/download_data.py      # UniProt 2026_03 snapshot   (~1 min)",
              "# Optional live refresh; overwrites the bundled 2026_03 snapshot (~1 min).\n"
              "# Skip this line to reproduce the paper from the bundled tables.\n"
              "# python src/download_data.py")
w = w.replace("`palm-dl/.cache/esm2_t12_35M_UR50D.pt` (the path `src/esm_embed.py` expects).",
              "`<workspace>/.cache/esm2_t12_35M_UR50D.pt`, i.e. a `.cache` directory next to the `palm-dl/`\n"
              "checkout (the path both `src/esm_embed.py` and `src/scan_embed.py` expect).")
wp.write_text(w, encoding="utf-8")

# 5) review package README/caveat template: author + candidate tiers wording
mp = ROOT / "src" / "make_review_package.py"
m = mp.read_text(encoding="utf-8")
m = m.replace("未并入最终发布模型", "未并入最终发布模型（序列-only 服务）")
m = m.replace("5% FDR 的 301 个候选位点", "5% 档的 301 个候选位点")
mp.write_text(m, encoding="utf-8")
print("eighth-round patch applied")
