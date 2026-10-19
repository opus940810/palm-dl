"""Patch remaining README sections and fix the subset count in the manuscript."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
readme = ROOT / "README.md"
lines = readme.read_text(encoding="utf-8").split("\n")


def replace_between(start_prefix: str, end_prefix: str, new: list[str]) -> None:
    global lines
    i = next(i for i, ln in enumerate(lines) if ln.startswith(start_prefix))
    j = next(j for j in range(i + 1, len(lines)) if lines[j].startswith(end_prefix))
    lines[i:j] = new


replace_between("AlphaFold pLDDT and C-alpha packing were extracted", "## Literature audit", [
    "AlphaFold pLDDT and C-alpha packing were extracted for the human proteins carrying annotated",
    "sites (1,950 sites, 463 positives in 476 proteins) and evaluated on the same cluster-disjoint",
    "folds:",
    "",
    "| Variant | AUC | delta |",
    "|---|---|---|",
    "| Forest, structured subset | 0.890 | — |",
    "| Forest + pLDDT + packing | 0.911 | +0.021 (4/5 folds) |",
    "",
    "The effect is biologically coherent — palmitoylated cysteines have lower pLDDT (66.2 vs 81.1)",
    "and lower packing (4.9 vs 10.8 C-alpha neighbours) — but modest, and it does not rescue the",
    "literature domain shift. The released model and service are therefore sequence-only; the paired",
    "analysis is reported in manuscript §2.12.",
    "",
])

readme.write_text("\n".join(lines), encoding="utf-8")

t = readme.read_text(encoding="utf-8")
t = t.replace("6,805 positive sites (477 experimental, 6,328 inferred), 10,567 hard negatives, 54,222\n"
              "background negatives, 23,897 proteins, 1,425 organisms, 64,143 unique 31-aa windows.",
              "4,424 positive sites (449 experimental, 3,975 inferred), 8,620 hard negatives, 54,222\n"
              "background negatives, 21,568 proteins, 786 organisms, 61,183 unique 31-aa windows.")
t = t.replace("PWM AUC\n0.537) and PalmDL separates it from same-protein cysteines only weakly (AUC 0.588; 0.607 against\nbackground)",
              "PWM AUC\n0.540) and PalmDL separates it from same-protein cysteines only weakly (AUC 0.585; 0.578 against\nbackground)")
readme.write_text(t, encoding="utf-8")

ms = ROOT / "manuscript" / "manuscript.md"
m = ms.read_text(encoding="utf-8")
m = m.replace("among all 26 non-empty member subsets by pooled", "among all 31 non-empty member subsets by pooled")
ms.write_text(m, encoding="utf-8")
print("readme/manuscript patched")
