"""Rewrite the manuscript for the S-only rebuild (v2).

Paragraphs are replaced wholesale by matching their first line; tables are
replaced by block between the table caption and the next blank line.
"""
from __future__ import annotations

import re
from pathlib import Path

MD = Path(__file__).resolve().parents[1] / "manuscript" / "manuscript.md"
t = MD.read_text(encoding="utf-8")
lines = t.split("\n")


def replace_paragraph(prefix: str, new: str) -> bool:
    global lines
    for i, ln in enumerate(lines):
        if ln.startswith(prefix):
            lines[i] = new
            return True
    print("NOT FOUND paragraph:", prefix[:60])
    return False


def replace_block(start_prefix: str, new_block: list[str]) -> bool:
    global lines
    for i, ln in enumerate(lines):
        if ln.startswith(start_prefix):
            j = i
            while j < len(lines) and lines[j].strip() != "":
                j += 1
            lines[i:j] = new_block
            return True
    print("NOT FOUND block:", start_prefix[:60])
    return False


replace_paragraph("**Results.**", (
    "**Results.** We assembled a versioned, evidence-stratified dataset from UniProt release 2026_03 "
    "comprising 67,266 labelled cysteine sites (4,424 S-palmitoyl cysteine sites, of which 449 are "
    "experimentally supported; 786 organisms; 557 human sites in 312 proteins) with explicitly "
    "separated hard (same-protein) and background negatives. Restricting positives to S-palmitoyl "
    "cysteine excluded 2,381 N-/O-palmitoyl annotations that a permissive text filter had admitted, "
    "and the corrected folds share no identical 31-residue window between training and test. Using a "
    "single model across four protocols we show that conventional benchmarks leak severely: 41.1% of "
    "test positives share an identical window with a training positive under a random site-level "
    "split and 45.2% under a protein-level split, inflating AUC from 0.944 to 0.985 and MCC from 0.720 "
    "to 0.810. PalmDL, a two-view ensemble of ESM-2 embeddings and motif/physicochemical features, "
    "reaches AUC 0.944 (95% CI 0.940–0.948), AUPR 0.730 and MCC 0.611 under cluster-held-out "
    "cross-validation; on held-out human sites it reaches 0.914 (0.900–0.927), marginally above the "
    "best published human value of 0.900, and 0.885 on a strict species holdout. It outperforms a "
    "CSS-Palm-style PWM, logistic regression, a multi-kernel CNN and both single views (ΔAUC +0.013). "
    "Probabilities are Platt-calibrated (ECE 0.025 at natural prevalence). A six-species transfer "
    "matrix shows that the code is shared almost losslessly across mammals and partly with plants, "
    "whereas bacterial lipidation is a different modification and is excluded by construction. "
    "In-silico mutagenesis recovers the cysteine-cluster rule, and FDR-controlled operating points "
    "support shortlists for experimental follow-up. A literature-mined external set of 136 sites from "
    "91 papers remains motif-atypical (AUC ≈0.59), which we report as a domain-shift limitation rather "
    "than as validation."
))

replace_paragraph("We queried UniProt release 2026_03", (
    "We queried UniProt release 2026_03 for all Swiss-Prot entries carrying at least one lipidation "
    "feature (9,089 entries) and kept only features whose note is **S-palmitoyl cysteine** on a "
    "cysteine residue. The first version of this filter matched any note containing \"palmitoyl\" and "
    "therefore admitted N-palmitoyl cysteine, N6-palmitoyl lysine and O-palmitoyl Ser/Thr; the "
    "corrected dataset contains **4,424 positive sites**, 449 supported by experimental evidence codes "
    "and 3,975 inferred (mostly by similarity, 89.9%), from 786 organisms. Human (557 sites in 312 "
    "proteins) and mouse (511) are the largest contributors (Figure 1a,b). Of the annotated positives, "
    "2,445 are cytoplasmic and only 27 extracellular (Figure 1c)."
))

replace_paragraph("Negatives were defined in two tiers.", (
    "Negatives were defined in two tiers. **Hard negatives** are unannotated cysteines of the same "
    "proteins whose ±15-residue windows do not overlap any palmitoylation annotation, whatever its "
    "chemistry (8,620 sites). **Background negatives** are cysteines sampled from 19,622 reviewed "
    "human proteins that carry no lipidation feature at all (54,222 sites). The final labelled set "
    "contains 67,266 cysteine sites from 21,568 proteins and 61,183 distinct 31-residue windows. "
    "Metrics remain high on the harder same-protein negatives alone (AUC 0.956; Table 3), so model "
    "ranking is not an artefact of easy background negatives."
))

replace_paragraph("We trained the same random forest on the same features under four protocols", (
    "We trained the same random forest on the same features under four protocols and evaluated each on "
    "holdout sites that were never used for fitting or threshold selection (Table 2, Figure 2). Under a "
    "random site-level split, **41.1% of test positives had an identical 31-residue window among the "
    "training positives**, and AUC reached 0.980 (MCC 0.843). A protein-level split, the most common "
    "precaution in the field, still leaked **45.2%** of test windows and gave AUC 0.985 (MCC 0.810). "
    "Cluster-level splitting removed the overlap entirely — **0% of any test row, positives and "
    "background negatives alike, shares a window with training** — and reduced AUC to 0.944 and MCC to "
    "0.720. Holding out an entire species gave AUC 0.885."
))

replace_paragraph("This audit has two consequences.", (
    "This audit has two consequences. First, headline numbers from site-level or protein-level splits "
    "cannot be compared with numbers from homology-controlled splits; a model reporting AUC > 0.97 is "
    "likely measuring peptide memorisation. Second, exact-window deduplication must cover the negative "
    "pool as well: in the first version of our splits, background proteins were not clustered and 751 "
    "windows (2,012 rows, almost all background negatives) crossed folds. Clustering every protein in "
    "the dataset and merging groups that share a window removes that residue. We therefore report every "
    "result in this paper under cluster-level splits, and additionally under a species holdout."
))

replace_paragraph("We compared five model families", (
    "We compared five model families under identical cluster-held-out folds (Table 1, Figure 3). A "
    "CSS-Palm-style position-specific log-odds matrix reached AUC 0.799 and MCC 0.344. Logistic "
    "regression on the flattened one-hot window reached AUC 0.829. A multi-kernel CNN reached AUC 0.854 "
    "and MCC 0.439. A random forest on one-hot and physicochemical channels reached AUC 0.931 and MCC "
    "0.595. An MLP on frozen ESM-2 (35M) embeddings reached AUC 0.931 and MCC 0.612. After the label "
    "correction the two strongest views are therefore tied on discrimination, with ESM-2 better at "
    "fixing a decision threshold."
))

replace_paragraph("The two strongest views make different errors", (
    "The two views make different errors: the forest captures exact positional chemistry and local motif "
    "syntax, whereas ESM-2 embeddings capture evolutionary and structural context. A rank-average of the "
    "two, followed by Platt calibration, was selected among all 26 non-empty member subsets by pooled "
    "validation AUC — never by test performance — and reached **AUC 0.944 (95% CI 0.940–0.948), AUPR "
    "0.730, MCC 0.611 and balanced MCC 0.766** (Table 1). Paired bootstrap over sites gives ΔAUC = "
    "+0.013 (95% CI 0.012–0.016) over either single view; the interval excludes zero. Adding the CNN to "
    "the ensemble lowered validation AUC (0.9420 versus 0.9487), so the two-view model was retained."
))

replace_paragraph("Because most published predictors are evaluated on human sites", (
    "Because most published predictors are evaluated on human sites, we repeated the analysis on "
    "human-only test folds (Table 3, Figure 4). With multi-species training and human test folds, PalmDL "
    "reaches **AUC 0.914 (0.900–0.927), AUPR 0.318, MCC 0.314 and balanced MCC 0.689**. A 4,000-"
    "replicate within-fold bootstrap against the best published human value, 0.900, gives P(AUC ≤ "
    "0.900) = 0.025: the estimate is marginally above the published value, but its interval touches it, "
    "so we describe PalmDL as matching, not surpassing, the state of the art. Training on human data "
    "alone is, by validation AUC, best served by the ESM-2 view alone (0.892; random forest 0.879), "
    "showing that non-human sites add little once homology is controlled. On the strict species "
    "holdout, in which no human protein and no cluster containing a human protein is seen during "
    "training and every human-propagated annotation is removed, the random forest reaches AUC 0.885 and "
    "balanced MCC 0.614. The 0.029 AUC gap between within-species and cross-species evaluation is much "
    "smaller than the 0.13 gap we measured before the label correction, because the N-palmitoyl "
    "bacterial sites that previously dominated the non-human training set are gone; what remains is a "
    "modest, measurable cost of species-specific context."
))

replace_paragraph("Probabilities from the calibrated ensemble are reliable", (
    "Probabilities from the calibrated ensemble are reliable at natural prevalence (ECE 0.025; Figure "
    "5), whereas the uncalibrated forest has ECE ≈ 0.12 and over-confident tails. On the hard-negative "
    "benchmark (positives plus same-protein cysteines only, 1:1.95 ratio), PalmDL reaches AUC 0.956 and "
    "MCC 0.678, compared with 0.810/0.414 for the PWM baseline (Table 3), so the advantage is not an "
    "artefact of easy background negatives."
))

replace_paragraph("Table 4 places our results next", (
    "Table 4 places our results next to the numbers reported by previous methods. Three patterns are "
    "visible. First, classical MCC values (0.45–0.67) were obtained under protocols that either reuse "
    "the training literature or do not hold out proteins, and the only classical tool with an "
    "independent test (CSS-Palm 2.0) dropped from MCC 0.671 to 0.523 when it did so. Second, the deep-"
    "learning systems with published AUCs (GPS-Palm 0.900, pCysMod 0.807) report neither cross-species "
    "nor explicit protein-disjoint evaluation. Third, the highest published number (Deep-Palm 0.931) "
    "comes from a site-level split, precisely the protocol that our audit shows to leak 41–45% of "
    "windows. PalmDL's 0.944 is obtained with zero shared windows, and its human estimate, 0.914, is "
    "marginally above the best published human AUC while being the only one accompanied by a leakage "
    "audit, a calibration analysis and a cross-species holdout. Because we could not rerun GPS-Palm or "
    "the retired classical tools on our folds, Table 4 compares protocols as much as algorithms; the "
    "released fold files remove this obstacle for future work."
))

replace_block("| **PalmDL** | this work | cluster-held-out CV, 6,805 sites", [
    "| **PalmDL** | this work | cluster-held-out CV, 4,424 sites | **AUC 0.944, MCC 0.611** | yes (0% shared windows) |",
    "| **PalmDL** | this work | human test folds | **AUC 0.914, balanced MCC 0.689** | yes (clusters) |",
    "| **PalmDL** | this work | species holdout | AUC 0.885, balanced MCC 0.614 | yes (strictest) |",
])

replace_paragraph("The species holdout in §2.4 removes every human-homologous cluster", (
    "The species holdout in §2.4 removes every human-homologous cluster and therefore measures transfer "
    "to genuinely novel families. To separate that from ordinary cross-species transfer we trained a "
    "separate PalmDL model on each of the six species with at least 85 S-palmitoyl cysteine sites "
    "(human, mouse, rat, cattle, pig, *A. thaliana*) and evaluated every model on every species, "
    "excluding target rows whose cluster appears in training (Figure 7, Table 5). Two regimes emerge. "
    "First, **mammals share one code**: models trained on human, mouse, rat, cattle or pig reach "
    "0.85–0.95 on the other mammals, comparable to the within-species diagonal. Second, **the code "
    "extends in part to plants**: mammalian models reach 0.96–0.98 on *A. thaliana* in the "
    "held-out-cluster setting, whereas the plant model reaches 0.70–0.85 on mammals — lower, but far "
    "above chance. Bacterial lipidation no longer enters the analysis, because after the label "
    "correction the bacterial positives are N-palmitoyl cysteines, a different modification that cannot "
    "be used to study the S-palmitoylation code."
))

replace_paragraph("First, **mammals share one code**.", (
    "The contrast with §2.4 is the key point. Cross-species transfer with orthologues available is "
    "essentially free (≈0.94), whereas transfer to human families absent from training is 0.885; most of "
    "the apparent generality of palmitoylation predictors therefore comes from orthologue coverage, not "
    "from an abstract motif grammar."
))

replace_block("| train \\ test | Human | Mouse | Rat | Cattle | Pig | *A. thaliana* | *E. coli* | *B. subtilis* |", [
    "| train \\ test | Human | Mouse | Rat | Cattle | Pig | *A. thaliana* |",
    "|---|---|---|---|---|---|---|",
    "| Human | {{M_HUMAN_HUMAN}} | {{M_HUMAN_MOUSE}} | {{M_HUMAN_RAT}} | {{M_HUMAN_CATTLE}} | {{M_HUMAN_PIG}} | {{M_HUMAN_ARAB}} |",
    "| Mouse | {{M_MOUSE_HUMAN}} | {{M_MOUSE_MOUSE}} | {{M_MOUSE_RAT}} | {{M_MOUSE_CATTLE}} | {{M_MOUSE_PIG}} | {{M_MOUSE_ARAB}} |",
    "| Rat | {{M_RAT_HUMAN}} | {{M_RAT_MOUSE}} | {{M_RAT_RAT}} | {{M_RAT_CATTLE}} | {{M_RAT_PIG}} | {{M_RAT_ARAB}} |",
    "| Cattle | {{M_CATTLE_HUMAN}} | {{M_CATTLE_MOUSE}} | {{M_CATTLE_RAT}} | {{M_CATTLE_CATTLE}} | {{M_CATTLE_PIG}} | {{M_CATTLE_ARAB}} |",
    "| Pig | {{M_PIG_HUMAN}} | {{M_PIG_MOUSE}} | {{M_PIG_RAT}} | {{M_PIG_CATTLE}} | {{M_PIG_PIG}} | {{M_PIG_ARAB}} |",
    "| *A. thaliana* | {{M_ARAB_HUMAN}} | {{M_ARAB_MOUSE}} | {{M_ARAB_RAT}} | {{M_ARAB_CATTLE}} | {{M_ARAB_PIG}} | {{M_ARAB_ARAB}} |",
    "",
    "Diagonal entries are from a protein-cluster-disjoint 20% holdout of the same species; off-diagonal entries use a model trained on all source rows, with target rows in the training clusters removed. Pig and *A. thaliana* have few held-out positives (19 and 18), so their diagonal is correspondingly noisy.",
])

MD.write_text("\n".join(lines), encoding="utf-8")
print("manuscript body updated; placeholders:",
      sum(1 for ln in lines if "{{M_" in ln))
