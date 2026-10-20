"""Second-round revisions: nested selection, cluster CIs, FDR wording, matrix v2."""
from __future__ import annotations

from pathlib import Path

MD = Path(__file__).resolve().parents[1] / "manuscript" / "manuscript.md"
lines = MD.read_text(encoding="utf-8").split("\n")


def rep(prefix: str, new: str) -> None:
    for i, ln in enumerate(lines):
        if ln.startswith(prefix):
            lines[i] = new
            return
    print("MISS:", prefix[:70])


def block(start: str, new: list[str]) -> None:
    for i, ln in enumerate(lines):
        if ln.startswith(start):
            j = i
            while j < len(lines) and lines[j].strip() != "":
                j += 1
            lines[i:j] = new
            return
    print("MISS block:", start[:70])


rep("**Results.** We assembled", (
    "**Results.** We assembled a versioned, evidence-stratified dataset from UniProt release 2026_03 "
    "comprising 67,266 labelled cysteine sites (4,424 S-palmitoyl cysteine sites, of which 455 are "
    "experimentally supported; 786 organisms; 557 human sites in 312 proteins) with explicitly "
    "separated hard (same-protein) and background negatives. Restricting positives to S-palmitoyl "
    "cysteine excluded 2,381 N-/O-palmitoyl annotations that a permissive text filter had admitted, and "
    "the corrected folds share no identical 31-residue window between training and test. Using a single "
    "model across four protocols we show that conventional benchmarks leak severely: 41.1% of test "
    "positives share an identical window with a training positive under a random site-level split and "
    "45.2% under a protein-level split, inflating AUC from 0.944 to 0.985 and MCC from 0.720 to 0.810. "
    "PalmDL, a two-view ensemble of ESM-2 embeddings and motif/physicochemical features whose member "
    "combination is chosen inside each outer fold, reaches AUC 0.938 (cluster-bootstrap 95% CI "
    "0.917–0.954), AUPR 0.697 and MCC 0.590 under cluster-held-out cross-validation; on held-out human "
    "sites it reaches 0.910 (0.890–0.929) and 0.885 on a strict species holdout. It improves on each "
    "single view (ΔAUC +0.008 over ESM-2, +0.007 over the forest, whose interval includes zero). "
    "Probabilities are Platt-calibrated (mean fold ECE 0.025). A six-species transfer matrix shows "
    "heterogeneous transfer rather than a universal code: most mammalian pairs transfer at 0.83–0.95, "
    "the plant model reaches 0.74–0.78 on mammals, and bacterial lipidation is excluded as a different "
    "modification. In-silico mutagenesis recovers the cysteine-cluster rule, and empirical "
    "validation-target operating points support candidate shortlists. A literature-mined external set "
    "of 136 sites from 91 papers remains motif-atypical (within-protein AUC 0.585), which we report as "
    "a domain-shift limitation rather than as validation."
))

rep("The two views make different errors", (
    "The two views make different errors: the forest captures exact positional chemistry and local motif "
    "syntax, whereas ESM-2 embeddings capture evolutionary and structural context. Member selection is "
    "now **nested inside each outer fold**: the 31 non-empty subsets are ranked by that fold's validation "
    "AUC, the winner is calibrated on the same validation rows, and only then applied to that fold's test "
    "rows. Four folds select RF+ESM and one selects LR+RF+CNN+ESM. The resulting ensemble reaches **AUC "
    "0.938 (cluster-bootstrap 95% CI 0.917–0.954), AUPR 0.697, MCC 0.590 and balanced MCC 0.747** (Table "
    "1). Paired cluster bootstrap gives ΔAUC = +0.008 (0.002–0.014) over ESM-2 and +0.007 "
    "(−0.004–0.018) over the forest, so the advantage over the forest is not resolved by these folds. For "
    "comparison, the previous globally selected combination scores 0.944 (0.926–0.957 cluster), which "
    "shows that pooling validation folds to choose one combination is mildly optimistic; we report the "
    "nested estimate as primary."
))

rep("Because most published predictors are evaluated on human sites", (
    "Because most published predictors are evaluated on human sites, we repeated the analysis on human-only "
    "test folds (Table 3, Figure 4). With multi-species training and human test folds, nested PalmDL "
    "reaches **AUC 0.910 (0.890–0.929), AUPR 0.294, MCC 0.297 and balanced MCC 0.662**. These values "
    "reproduce the published human AUC of 0.900 only descriptively: no paired test is possible across "
    "different test sets, and the within-study bootstrap probability that our estimate falls below 0.900 "
    "is 0.08. We therefore report the number and leave the comparison to context. Training on human data "
    "alone is selected as RF+ESM in four folds and ESM alone in one (mean 0.904). On the strict species "
    "holdout the random forest reaches AUC 0.885 and balanced MCC 0.614. The 0.025 difference between "
    "the human CV estimate and the species holdout uses different model families and negative "
    "compositions, so it is descriptive rather than a clean decomposition of a species-specific cost."
))

rep("Probabilities from the calibrated ensemble are reliable", (
    "Probabilities from the calibrated ensemble are reliable at natural prevalence (mean fold ECE 0.025, "
    "pooled 10-bin ECE 0.021; Figure 5), whereas the uncalibrated forest has ECE ≈ 0.12 and over-confident "
    "tails. On the hard-negative benchmark (positives plus same-protein cysteines only, 1:1.95 ratio), "
    "PalmDL reaches AUC 0.952 and MCC 0.664, compared with 0.811/0.414 for the PWM baseline (Table 3), so "
    "the advantage is not an artefact of easy background negatives."
))

rep("Table 4 places our results next", (
    "Table 4 places our results next to the numbers reported by previous methods. Three patterns are "
    "visible. First, classical MCC values (0.45–0.67) were obtained under protocols that either reuse the "
    "training literature or do not hold out proteins, and the only classical tool with an independent test "
    "(CSS-Palm 2.0) dropped from MCC 0.671 to 0.523 when it did so. Second, the deep-learning systems with "
    "published AUCs (GPS-Palm 0.900, pCysMod 0.807) report neither cross-species nor explicit "
    "protein-disjoint evaluation. Third, the highest published number (Deep-Palm 0.931) comes from a "
    "site-level split, precisely the protocol that our audit shows to leak 41–45% of windows. PalmDL's "
    "0.938 is obtained with zero shared windows under nested member selection, and its human estimate, "
    "0.910, is reported for context rather than as a paired superiority claim. Because we could not rerun "
    "GPS-Palm or the retired classical tools on our folds, Table 4 compares protocols as much as "
    "algorithms; the released fold files remove this obstacle for future work."
))

block("| **PalmDL** | this work | cluster-held-out CV, 4,424 sites", [
    "| **PalmDL** | this work | cluster-held-out CV, 4,424 sites | **AUC 0.938, MCC 0.590** | yes (0% shared windows) |",
    "| **PalmDL** | this work | human test folds | **AUC 0.910, balanced MCC 0.662** | yes (clusters) |",
    "| **PalmDL** | this work | species holdout | AUC 0.885, balanced MCC 0.614 | yes (strictest) |",
])

block("| train \\ test | Human | Mouse | Rat | Cattle | Pig | A. thaliana |", [
    "| train \\ test | Human | Mouse | Rat | Cattle | *A. thaliana* | Pig |",
    "|---|---|---|---|---|---|---|",
    "| Human | 0.917 | 0.928 | 0.910 | 0.948 | 0.980 | — |",
    "| Mouse | 0.899 | 0.877 | 0.880 | 0.883 | 0.984 | — |",
    "| Rat | 0.890 | 0.914 | 0.928 | 0.945 | 0.984 | 0.940 |",
    "| Cattle | 0.899 | 0.903 | 0.917 | 0.880 | 0.979 | 0.937 |",
    "| *A. thaliana* | 0.739 | 0.766 | 0.781 | 0.775 | — | 0.785 |",
    "| Pig | 0.834 | 0.854 | 0.858 | 0.880 | 0.964 | — |",
    "",
    "One model per source species is trained on its positives and same-protein negatives (no human background); 20% of source clusters are held out for the diagonal and another 15% for early stopping. Target cells use only positives and same-protein negatives, and rows whose cluster appears in either the fit or the validation split are removed. Off-diagonal cells therefore share the model of their row, as stated in Methods. Dashes mark cells with fewer than 20 eligible test positives (Human→Pig 13, Mouse→Pig 19, and the *A. thaliana* and Pig diagonals with 15 and 16).",
])

rep("The species holdout in §2.4 removes every human-homologous cluster", (
    "The species holdout in §2.4 removes every human-homologous cluster and therefore measures transfer to "
    "genuinely novel families. To separate that from ordinary cross-species transfer we trained one model "
    "per species on its positives and same-protein negatives, holding out clusters for the diagonal and for "
    "early stopping, and evaluated every model on every species with target rows restricted to positives "
    "and same-protein negatives (Figure 7, Table 5). Transfer is **heterogeneous**, not uniformly "
    "lossless. Most mammalian pairs fall between 0.83 and 0.95 (mean 0.896), with Pig→Human (0.834) the "
    "clear exception, 0.11 below the human diagonal and 0.14 below the pig diagonal. The plant model "
    "reaches 0.74–0.78 on mammals, whereas mammalian models reach 0.96–0.98 on *A. thaliana*. Several "
    "cells have too few eligible positives after cluster removal and are left empty rather than "
    "estimated from a handful of sites. The matrix therefore describes transfer among these sampled "
    "species and negative-sampling protocols; it does not isolate the causal contribution of orthologue "
    "coverage."
))

rep("The contrast with §2.4 is the key point.", (
    "Two consequences follow. First, claims of a single shared eukaryotic code are not supported by this "
    "design: the pig exception and the plant asymmetry are as informative as the mammalian mean. Second, "
    "the strict human holdout (0.885) and the human CV estimate (0.910) differ by only 0.025 and use "
    "different model families and negative compositions, so they should not be read as a clean "
    "species-specific penalty. A paired design with matched negatives and a frozen target set would be "
    "needed to separate orthologue coverage from motif generalisation."
))

rep("Because a predictor is only useful", (
    "Because a predictor is only useful if its errors are predictable, we stratified the held-out predictions "
    "(Table 6). Discrimination remains better for inferred than for experimentally supported sites (AUC "
    "0.947 versus 0.907). Evidence tiers now follow one rule: a site counts as experimental when **any** "
    "evidence clause is experimental (ECO:0000269, 0007744, 0007005 or 0006056), which gives 455 "
    "experimental and 3,969 inferred sites; six mixed-clause sites were previously misclassified by a "
    "comma-only split and are now aligned with Table 6. Performance is stable across annotated topologies "
    "(cytoplasmic 0.887, extracellular 0.855, topology unknown 0.898), with the lumenal class too small to "
    "interpret (14 positives). Topology is itself an input feature, so these strata describe where the "
    "model is accurate, not an independent validation of the topology signal."
))

rep("For practical use we derived", (
    "For practical use we derived, **within each fold**, the lowest threshold whose complete score groups "
    "satisfy an empirical validation false-discovery target, and evaluated it on that fold's outer test "
    "rows (Figure 8b). The nominal targets are not achieved out of sample, and we report both sides: at a "
    "1% validation target the held-out precision is 0.953, recall 0.201 and observed false-discovery "
    "proportion 4.7%; at 5%, 0.900/0.331/10.0%; at 10%, 0.849/0.398/15.1%; at 20%, 0.759/0.490/24.1%. "
    "These are empirically calibrated operating points, not a demonstrated FDR-controlling procedure; the "
    "gap between validation target and held-out FDP is the honest measure of their transferability."
))

rep("It is not a clean benchmark.", (
    "It is not a clean benchmark. The motif that the model learned is barely enriched in this set (PWM AUC "
    "0.537 on the 133 scorable sites), and PalmDL separates the literature sites from same-protein "
    "cysteines only weakly (weighted within-protein AUC 0.585; pooled novel-site AUC 0.578). Manual "
    "inspection explains why: genuine sites are mixed with extraction errors, and many genuine sites are "
    "transmembrane-proximal or otherwise motif-atypical cysteines that an annotation-trained model has "
    "least often seen. A low PWM score does not by itself prove that a literature site is wrong; it shows "
    "that this class lies outside the motif distribution the model was trained on."
))

rep("To make sure this was not an extraction artefact", (
    "To make sure this was not an extraction artefact, we curated all 136 sites by reading every sentence. "
    "Ninety-nine mentions could be assigned unambiguously to a protein and a cysteine — either because the "
    "sentence states the site explicitly or because a cysteine-to-serine/alanine mutant lost palmitoylation "
    "— while 37 were rejected as sites of a different protein in the same sentence, catalytic-site mutants "
    "of DHHC enzymes, predicted-only sites, negative results, or table rows with mismatched proteins. The "
    "curated subset behaves like the full set, so extraction noise is not the explanation."
))

rep("Consistent with this, the curated sites", (
    "The curated sites are almost absent from the candidate list: **none of the 99 appears at the 5% tier, "
    "one at the 20% tier and two at the 30% tier** (the same counts as the full 136-site set). The honest "
    "conclusion is a domain-shift limitation that applies to any annotation-trained predictor, including "
    "PalmDL: it reproduces and extends the UniProt consensus, but does not recognise the motif-atypical, "
    "transmembrane-proximal sites that mechanistic papers tend to characterise."
))

block("| All verified sites | 136 | 0.540 | 0.585 |", [
    "| All verified sites | 133 | 0.537 | 0.588 |",
])

rep("The 5% tier is biologically coherent.", (
    "The 5% tier is biologically coherent in composition. 61.8% of its sites carry annotated cytoplasmic "
    "topology versus 5.6% of all scanned cysteines (11.1-fold). The highest-confidence examples at the 5% "
    "tier include endothelin receptor A C386 (q = 0.044); its neighbour C388 is in the 20% tier "
    "(q = 0.174), and the other candidates mentioned in earlier drafts (claudin-2 C185, mucolipin-3 C550) "
    "have q = 0.70 and 0.41 and are not top-tier examples. Topology is an input feature, so this "
    "enrichment describes the composition of the list rather than validating the predictions."
))

rep("Three caveats bound the resource.", (
    "Four caveats bound the resource. First, the proteome-wide prior is small (557 annotated sites among "
    "262,286 cysteines, ≈0.2%), so any strict operating point is intrinsically severe. Second, 82% of "
    "scanned cysteines have no topology annotation and are strongly depleted from the shortlist (38.2%), "
    "so the list is biased toward proteins with curated topology. Third, the decoy tail is thin — few "
    "shuffled decoys reach the top of the score distribution — so the FDR estimate is optimistic and the "
    "tiers are candidate rankings rather than proven error rates. Fourth, known-site recovery in the scan "
    "(44.5% at the 5% tier) is an in-sample upper bound because the released model was trained on all "
    "UniProt annotations; the corresponding out-of-fold recall is 33.1% (§2.8)."
))

rep("On those folds PalmDL reaches AUC 0.944 overall", (
    "On those folds nested PalmDL reaches AUC 0.938 overall and 0.910 on human sites. The best previously "
    "reported human AUC is 0.900 (GPS-Palm [8]); the best overall figure in the literature is 0.931 from a "
    "preprint using a site-level split (Deep-Palm [14]). These are different test sets and negative "
    "schemes, so the numbers are context rather than a paired comparison; the within-study bootstrap "
    "probability that our human estimate falls below 0.900 is 0.08. This is not a head-to-head rerun: "
    "GPS-Palm's web server has been retired (HTTP 404 on 18 October 2026) and its standalone distribution "
    "is a GUI-only academic-licence package. A true head-to-head on our released folds remains the single "
    "most important piece of missing evidence, and we provide the fold files to make it straightforward."
))

rep("The cross-species experiments split into two statements", (
    "The cross-species experiments are deliberately descriptive. With orthologues available, most "
    "mammalian pairs transfer at 0.83–0.95, but Pig→Human is 0.834 and the plant→mammal direction is "
    "0.74–0.78; without human-homologous clusters, human performance is 0.885. The matrix therefore shows "
    "heterogeneous transfer among sampled species and negative-sampling protocols, not a uniform "
    "eukaryotic code, and it does not isolate orthologue coverage as the cause. Bacterial lipidation is "
    "excluded rather than compared: after the label correction the bacterial positives are N-palmitoyl "
    "cysteines, a different modification."
))

rep("Two design choices follow from the error analysis.", (
    "Two design choices follow from the error analysis. First, the gap between inferred (AUC 0.947) and "
    "experimental sites (0.907) is an annotation-memory effect that any user should know about. Second, "
    "reporting a fixed 0.5 cut-off is less useful than reporting empirical operating points, provided "
    "their transfer is shown: at a 5% validation target the held-out precision is 0.900 and the observed "
    "false-discovery proportion 10.0%, so the nominal target is not a guarantee."
))

rep("Structure is a modest, real improvement rather than a cure.", (
    "Structure is a modest, real improvement rather than a cure. Three AlphaFold-derived numbers — pLDDT "
    "at the site, local pLDDT and C-alpha packing — improve the forest view from 0.890 to 0.911, "
    "consistently in four of five folds. pLDDT is a model-confidence measure that correlates with "
    "disorder, not a direct experimental measure of flexibility. The domain shift documented in §2.10 is "
    "not simply a burial-versus-exposure problem, and structure alone will not close it."
))

# limitations
rep("**Limitations.** (i) 90% of positives", (
    "**Limitations.** (i) 90% of positives are inferred rather than experimental; we stratify by evidence "
    "but cannot eliminate label noise, and the 455 experimental sites are too few for a well-powered "
    "experimental-only benchmark. (ii) Background negatives are proteins without annotation, not proteins "
    "proven not to be palmitoylated. (iii) The cluster definition uses 4-mer cosine similarity, an "
    "operational proxy: cosine 0.28–0.34 corresponded to median 71% global identity in 30 sampled pairs, "
    "so remote homologues below ~50% identity may still cross folds. (iv) We model only cysteine-centred "
    "windows and do not use predicted structure, DHHC-enzyme identity, or tissue expression. (v) The "
    "ensemble is calibrated at natural prevalence; performance at other prevalences requires "
    "recalibration. (vi) Independent, site-level public benchmarks do not exist: SwissPalm is behind an "
    "IP-verification gate, dbPTM and CysModDB no longer resolve, the MDD-Palm supplements contain no site "
    "lists, and GPS-Palm's server is retired. Both alternatives we could obtain fail a motif sanity check, "
    "so neither is used as a benchmark. (vii) The proteome-wide scan uses the forest view and a "
    "shuffled-flank null with a thin tail; the tiers are candidate rankings, not demonstrated error "
    "rates. (viii) The cross-species matrix samples heterogeneous negative sets and leaves four cells "
    "empty for lack of eligible positives; it does not support a causal orthologue-coverage claim. (ix) "
    "The deployed service averages view logits rather than within-fold ranks, which on the same folds "
    "gives AUC 0.940 versus 0.944 for the rank average under the global combination."
))

# methods
rep("Lipidation features were parsed with a regular expression", (
    "Lipidation features were parsed with a regular expression. A site was labelled positive only when the "
    "note was **S-palmitoyl cysteine** and the canonical residue was a cysteine. A site was classified as "
    "experimental when at least one ECO clause in its evidence string matched ECO:0000269, ECO:0007744, "
    "ECO:0007005 or ECO:0006056; clauses are separated by ',' or ';' in UniProt, and both separators are "
    "handled. This gives 455 experimental and 3,969 inferred sites, the counts used in Table 6. Each site "
    "was represented by a 31-residue window centred on the cysteine and zero-padded with `X`. Hard "
    "negatives were unannotated cysteines of positive-carrying proteins at least 16 residues from any "
    "palmitoylation site, capped at four per protein; background negatives were up to three cysteines from "
    "each human protein without any lipidation feature. Topology was derived from UniProt topological "
    "domains."
))

rep("Every protein in the dataset", (
    "Every protein in the dataset — positive-carrying and background alike — was clustered by cosine "
    "similarity of its 4-mer count vector at 0.30, and clusters sharing an identical 31-residue window were "
    "merged, so no exact peptide occurs in two folds. Clusters were assigned to five folds with greedy "
    "balancing of positive counts (885 positives per fold). The species holdout used human positives and "
    "human same-protein negatives as test; background negatives are human by construction and were excluded "
    "from both sides. Training used non-human sites after removing every cluster containing a human protein "
    "and every annotation propagated from a human accession, leaving 1,480 training and 557 test positives."
))

rep("For each outer fold, 15% of the training clusters", (
    "For each outer fold, 15% of the training clusters were held out as inner validation. Background "
    "negatives were subsampled to 1× the number of positives in training; test folds retained their natural "
    "composition. Thresholds were chosen to maximise MCC on the inner validation set. Ensemble member "
    "selection, calibration and thresholding were performed **within the same outer fold**: the 31 non-empty "
    "subsets were ranked by that fold's validation AUC, the winner was Platt-calibrated on those rows, and "
    "the calibrated scores were applied to that fold's test rows. The globally selected combination is "
    "reported only as a sensitivity analysis. Discrimination is the unweighted mean of the five fold "
    "metrics; MCC and accuracy pool confusion counts at the per-fold thresholds. Two bootstrap schemes are "
    "reported: resampling sites and resampling protein clusters within folds."
))

rep("Six species with at least 85 S-palmitoyl cysteine sites", (
    "Six species with at least 85 S-palmitoyl cysteine sites were used: human, mouse, rat, cattle, pig and "
    "*Arabidopsis thaliana*. For each source species one model (random forest plus ESM head) was trained on "
    "that species' positives and its same-protein negatives; human background negatives were excluded. "
    "Twenty percent of the source clusters were held out for the diagonal and another 15% for early "
    "stopping. Off-diagonal cells use the same model as their row. Target rows were restricted to positives "
    "and same-protein negatives and removed if their cluster appeared in either the fit or the validation "
    "split. Cells with fewer than 20 eligible test positives are reported as empty rather than estimated."
))

rep("For each target FDR f and each fold separately", (
    "For each empirical target f and each fold separately, validation rows were sorted by predicted "
    "probability and the lowest threshold was chosen whose **complete score groups** satisfied the "
    "cumulative false-discovery constraint (FP / (TP + FP) ≤ f); splitting a tied group is not allowed. The "
    "threshold was applied to that fold's outer test rows and precision, recall, predicted counts and the "
    "observed test false-discovery proportion were pooled across folds. Validation targets are nominal: "
    "they are not a guarantee of out-of-sample FDR control, and both values are reported."
))

rep("All 262,286 cysteines of the reviewed human sequences", (
    "All 262,286 cysteines of the reviewed human sequences were scored in chunks with the forest view of the "
    "released model, calibrated on the same out-of-fold validation rows. The ESM view was not run at this "
    "scale. The null was 20,000 shuffled-flank decoys that keep the central cysteine but destroy the "
    "motif; a C→A decoy would remain motif-shaped. FDR at a threshold was estimated as (decoys at or above "
    "the threshold) × (targets / decoys) divided by (targets at or above the threshold), with q-values as "
    "the running minimum over lower thresholds. Because the decoy tail is thin and the shuffled null is not "
    "proven to represent unmodified cysteines, the tiers are candidate rankings rather than demonstrated "
    "error rates."
))

rep("Human proteins with at least one annotated site", (
    "Human proteins with at least one annotated site and a random sample of background human proteins were "
    "queried through the AlphaFold DB prediction API; models were retrieved for chains shorter than 700 "
    "residues, covering 1,950 human sites (463 positives in 476 proteins). Per site we recorded the C-alpha "
    "pLDDT, mean pLDDT within ±5 positions, and C-alpha atoms within 10 Å excluding neighbours closer than "
    "four positions. The same random forest and cluster-disjoint folds were used with and without the "
    "features, so the comparison is paired at the fold level. The prediction service exposes the sequence "
    "path only."
))

# figure legends
rep("**Figure 3. Discrimination under cluster-held-out cross-validation.**", (
    "**Figure 3. Discrimination under cluster-held-out cross-validation.** (a) ROC and (b) precision–recall "
    "curves pooled over out-of-fold predictions; the table reports unweighted fold means, so the pooled "
    "curves and the table are related but not identical estimands. Models: PWM, logistic regression, random "
    "forest, CNN, ESM-2 view and nested PalmDL (4,424 positives, 67,266 sites)."
))
rep("**Figure 4. Performance on human sites under three evaluation regimes.**", (
    "**Figure 4. Performance on human sites under three evaluation regimes.** (a) AUC and (b) balanced MCC "
    "for multi-species training with human test folds, human-only training, and the strict species holdout. "
    "The species-holdout bars are single-model (random forest/CNN) results and are not a full ensemble; "
    "negative compositions differ between the three regimes."
))
rep("**Figure 7. The palmitoylation code is shared across mammals and partly with plants.**", (
    "**Figure 7. Cross-species transfer is heterogeneous.** AUC of single-species PalmDL models, rows = "
    "training species, columns = test species, for the six species with at least 85 S-palmitoyl cysteine "
    "sites. Most mammalian pairs transfer at 0.83–0.95 with a pig→human exception; the plant model reaches "
    "0.74–0.78 on mammals. Empty cells have fewer than 20 eligible test positives after cluster removal. "
    "Diagonal entries come from a protein-cluster-disjoint 20% holdout of the same species."
))

# table 5 caption / table 8 caption wording
rep("**Table 5. Cross-species transfer matrix (AUC) for six species.**", "**Table 5. Cross-species transfer matrix (AUC) for six species.**")
rep("**Table 7. Automated literature mining as an external benchmark: yields and motif consistency.**", (
    "**Table 7. Automated literature mining as an external benchmark.**"
))

MD.write_text("\n".join(lines), encoding="utf-8")
print("second-round manuscript patch applied")
