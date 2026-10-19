"""Third pass: fix duplicated paragraphs, rewrite the scan section, update table 8."""
from __future__ import annotations

from pathlib import Path

MD = Path(__file__).resolve().parents[1] / "manuscript" / "manuscript.md"
lines = MD.read_text(encoding="utf-8").split("\n")

out = []
removed = 0
for ln in lines:
    if ln.startswith("The contrast with §2.4 is the key point. Cross-species transfer with orthologues available is essentially free (0.94)"):
        removed += 1
        continue
    if ln.startswith("Diagonal entries are from a protein-cluster-disjoint 20% holdout of the same species; off-diagonal entries use a model trained on all rows"):
        removed += 1
        continue
    out.append(ln)
lines = out
print("removed duplicates:", removed)


def rep(prefix: str, new: str) -> None:
    for i, ln in enumerate(lines):
        if ln.startswith(prefix):
            lines[i] = new
            return
    print("MISS:", prefix[:60])


def block(start: str, new: list[str]) -> None:
    for i, ln in enumerate(lines):
        if ln.startswith(start):
            j = i
            while j < len(lines) and lines[j].strip() != "":
                j += 1
            lines[i:j] = new
            return
    print("MISS block:", start[:60])


rep("### 2.7 The palmitoylation code is eukaryotic, not species-specific",
    "### 2.7 The palmitoylation code is shared across mammals and partly with plants")

rep("As a forward-looking counterpart to the retrospective benchmarks", (
    "As a forward-looking counterpart to the retrospective benchmarks, we applied the forest view of the "
    "released model to every cysteine of the reviewed human proteome — **262,286 sites in 19,814 "
    "proteins** — and controlled the false-discovery rate with 20,000 shuffled-flank decoys that keep the "
    "central cysteine but destroy the motif. (A C→A decoy is inappropriate: §2.9 shows that the model is "
    "deliberately insensitive to the identity of the centre.) At the 1% tier the scan yields 215 sites in "
    "134 proteins; at 5%, 301 sites in 199 proteins (248 already annotated, 53 novel); at 10%, 396; at "
    "20%, 615; and at 30%, 1,010 (Table 8, Figure 9)."
))

rep("The shortlist is biologically coherent.", (
    "The 5% tier is biologically coherent. 61.8% of its sites carry annotated cytoplasmic topology versus "
    "5.6% of all scanned cysteines — an 11.1-fold enrichment — and its highest-scoring entries include "
    "juxtamembrane cysteines of single-pass membrane proteins (endothelin receptor A C386/C388, claudin-2 "
    "C185, mucolipin-3 C550). Only one of the 136 literature-mined sites reaches the 20% tier, consistent "
    "with the motif-atypical nature of that set (§2.10)."
))

rep("Two caveats bound the resource.", (
    "Three caveats bound the resource. First, the proteome-wide prior is small (557 annotated sites among "
    "262,286 cysteines, ≈0.2%), so strict FDR control is intrinsically severe; users interested in "
    "discovery should work at the 20–30% tiers. Second, 82% of scanned cysteines have no topology "
    "annotation and are strongly depleted from the shortlist (38.2%), so the list is biased toward "
    "proteins with curated topology. Third, and most important, the decoy tail is thin — few shuffled "
    "decoys reach the top of the score distribution — so the estimated FDR is optimistic and the tiers "
    "should be read as candidate rankings rather than as proven error rates. Known-site recovery in the "
    "scan (44.5% at the 5% tier) is also an in-sample upper bound, because the released model was trained "
    "on all UniProt annotations; the corresponding out-of-fold recall is 34.5% (§2.8)."
))

block("| 1% | 217 | 138 | 71 | 12.4% | 146 | 0.969 |", [
    "| 1% | 215 | 134 | 197 | 35.4% | 18 | 0.979 |",
    "| 5% | 301 | 199 | 248 | 44.5% | 53 | 0.959 |",
    "| 10% | 396 | 269 | 276 | 49.6% | 120 | 0.938 |",
    "| 20% | 615 | 416 | 334 | 60.0% | 281 | 0.894 |",
    "| 30% | 1,010 | 701 | 380 | 68.2% | 630 | 0.830 |",
])

rep('"Known" refers to the 572 S-palmitoyl cysteine sites', (
    "\"Known\" refers to the 557 S-palmitoyl cysteine sites annotated in UniProt 2026_03 for human "
    "proteins; \"novel\" sites are unannotated cysteines. FDR is estimated by target–decoy comparison "
    "against 20,000 shuffled-flank decoys scored with the forest view; the decoy tail is thin, so these "
    "tiers are candidate rankings rather than proven error rates."
))

MD.write_text("\n".join(lines), encoding="utf-8")
print("placeholders left:", sum(1 for ln in lines if "{{" in ln))
