"""Second pass: error structure, label-bug retraction, literature, structure,
discussion, limitations, methods and legends."""
from __future__ import annotations

from pathlib import Path

MD = Path(__file__).resolve().parents[1] / "manuscript" / "manuscript.md"
lines = MD.read_text(encoding="utf-8").split("\n")


def rep(prefix: str, new: str) -> None:
    for i, ln in enumerate(lines):
        if ln.startswith(prefix):
            lines[i] = new
            return
    print("MISS paragraph:", prefix[:60])


def block(start: str, new: list[str]) -> None:
    for i, ln in enumerate(lines):
        if ln.startswith(start):
            j = i
            while j < len(lines) and lines[j].strip() != "":
                j += 1
            lines[i:j] = new
            return
    print("MISS block:", start[:60])


rep("Because a predictor is only useful", (
    "Because a predictor is only useful if its errors are predictable, we stratified the held-out "
    "predictions (Table 6). Discrimination remains better for inferred than for experimentally "
    "supported sites (AUC 0.947 versus 0.907), the annotation-memory effect expected when 90% of labels "
    "are propagated: experimentally characterised sites are the atypical, family-specific cases that "
    "propagation does not cover, and they are exactly the ones a user is most likely to care about. "
    "Performance is stable across annotated topologies (cytoplasmic 0.887, extracellular 0.855, "
    "topology unknown 0.898), with the lumenal class too small to interpret (14 positives)."
))

rep("For practical use we derived thresholds", (
    "For practical use we derived, **within each fold**, the threshold whose estimated false-discovery "
    "rate is at most a target value, and evaluated it on that fold's outer test rows (Figure 8b). At an "
    "estimated FDR of 1% the threshold retains 20.2% of known sites at 95.3% precision; at 5%, 34.5% at "
    "90.3%; at 10%, 40.4% at 85.5%; at 20%, 53.5% at 76.6%. Thresholds are deliberately conservative, "
    "because the proteome-wide prior is small, and they give a principled way to generate a shortlist "
    "without committing to an arbitrary score cut-off."
))

block("| Experimentally supported sites | 506", [
    "| Experimentally supported sites | 455 | 0.907 | 0.176 |",
    "| Inferred sites | 3,969 | 0.947 | 0.682 |",
    "| Cytoplasmic | 2,445 | 0.887 | 0.829 |",
    "| Extracellular | 27 | 0.855 | 0.375 |",
    "| Topology unknown | 1,801 | 0.898 | 0.501 |",
])

rep("In-silico mutagenesis of 300 human positive windows", (
    "In-silico mutagenesis of 300 human positive windows (every position mutated to every amino acid, "
    "scored with the forest view; Figure 8a) shows that substitutions at the central position are "
    "neutral, as they must be: every candidate is a cysteine, so the centre carries no discriminative "
    "information. The informative positions are the flanks. Introducing a cysteine at +1 increases the "
    "predicted probability by 0.089 and at −1 by 0.083, recovering the known preference for cysteine "
    "clusters; introducing glutamate or aspartate at −1 decreases it by 0.06–0.09, recovering the "
    "exclusion of acidic residues near the modification site."
))

rep("The same audit exposed a label-quality problem.", (
    "The same audit exposed a bug in **our own label filter**, not in UniProt. The first version of the "
    "dataset matched any lipidation note containing \"palmitoyl\", which admitted N-palmitoyl cysteine "
    "(2,315 sites), N6-palmitoyl lysine (34), N-palmitoyl glycine (18) and O-palmitoyl serine/threonine "
    "(14) — 35% of the original positives. These are legitimate annotations of different chemical "
    "modifications; they were mis-classified by us, and the earlier claim that they reflected UniProt "
    "annotation errors was wrong and is withdrawn. The corrected dataset requires an S-palmitoyl "
    "cysteine note on a cysteine residue and contains no non-cysteine positives. The label correction "
    "changed every downstream number in this paper; the audit trail (per-site notes before and after) is "
    "released with the data, and we report the corrected results throughout."
))

rep("It is not a clean benchmark.", (
    "It is not a clean benchmark. The motif that the model learned is barely enriched in this set (PWM "
    "AUC 0.540 against same-protein cysteines), and PalmDL separates the literature sites from those "
    "cysteines only weakly (AUC 0.585). Manual inspection explains why: genuine sites are mixed with "
    "extraction errors, and many genuine sites (MyD88 C113, iRhom2 C476, GRP78 C420) are transmembrane-"
    "proximal cysteines whose linear motif is weak, which is precisely the class that an annotation-"
    "trained model has least often seen."
))

rep("To make sure this was not an extraction artefact", (
    "To make sure this was not an extraction artefact, we curated all 136 sites by reading every "
    "sentence. Ninety-nine mentions could be assigned unambiguously to a protein and a cysteine — either "
    "because the sentence states the site explicitly or because a cysteine-to-serine/alanine mutant lost "
    "palmitoylation — while 37 were rejected as sites of a different protein in the same sentence, "
    "catalytic-site mutants of DHHC enzymes, predicted-only sites, negative results, or table rows with "
    "mismatched proteins. The curated subset behaves like the full set, so extraction noise is not the "
    "explanation. Allowing every site and every same-protein negative to score its best window within "
    "±2 residues also leaves the ranking unchanged, so numbering offsets do not explain the failure "
    "either."
))

rep("Consistent with this, the curated sites", (
    "Consistent with this, the curated sites are essentially absent from the proteome-wide shortlist: one "
    "of 99 appears at the 5% tier (expected 0.12 by chance) and two at the 20% tier. The honest "
    "conclusion is a domain-shift limitation that applies to any annotation-trained predictor, including "
    "PalmDL: it reproduces and extends the UniProt consensus, but does not recognise the motif-atypical, "
    "transmembrane-proximal sites that mechanistic papers tend to characterise."
))

block("| All verified sites | 133", [
    "| All verified sites | 136 | 0.540 | 0.585 |",
])

rep("Sequence-only models cannot see whether", (
    "Sequence-only models cannot see whether a cysteine is buried in a folded domain or exposed on a "
    "membrane-proximal surface, which matters because the acyl chain must be accommodated. We fetched "
    "AlphaFold models for the human proteins carrying annotated sites and computed three features per "
    "cysteine: pLDDT at the site, mean pLDDT within ±5 residues, and C-alpha packing density (number of "
    "C-alpha atoms within 10 Å, excluding sequence neighbours). 1,950 human sites (463 positives in 476 "
    "proteins) were covered."
))

rep("We compared four ways of using them", (
    "Adding the features to the random-forest view lifts AUC from 0.890 to 0.911 on the structured "
    "subset (Δ +0.021, positive in 4 of 5 folds; Table 9). The effect is modest because the language-"
    "model view already captures much of the structural context, and because only annotated human sites "
    "have models. We therefore report structure as a paired analysis rather than as part of the released "
    "sequence model."
))

rep("The direction of the structural effect is biologically coherent", (
    "The direction of the structural effect is biologically coherent: palmitoylated cysteines have lower "
    "pLDDT (66.2 versus 81.1) and lower packing (4.9 versus 10.8 C-alpha neighbours), i.e. they sit in "
    "less ordered, more exposed segments where the acyl chain can be accommodated and where DHHC enzymes "
    "act."
))

rep("It does not cure the domain shift", (
    "It does not cure the domain shift of §2.10. The class of sites that literature studies characterise "
    "— often transmembrane-proximal, sometimes weakly structured — is not rescued by burial or "
    "flexibility alone. Structure therefore improves ranking where annotations already exist, but it is "
    "not by itself the missing piece for genuinely novel sites."
))

block("| Forest + structure | 0.890", [
    "| Forest + structure | 0.890 → 0.911 (1,950 sites) | 0.911 | +0.021 (4/5 folds) |",
])

block("| Feature | Positive sites (mean) | Negative sites (mean) |", [
    "| Feature | Positive sites (mean) | Negative sites (mean) |",
    "|---|---|---|",
    "| pLDDT at cysteine | 66.2 | 81.1 |",
    "| C-alpha neighbours within 10 Å | 4.9 | 10.8 |",
])

rep("We set out to build a palmitoylation predictor", (
    "We set out to build a palmitoylation predictor that could credibly claim state-of-the-art "
    "performance. Doing so required first showing that the field's benchmarks overstate performance. With "
    "a single fixed model we measured that random and protein-level splits share 41–45% of test-positive "
    "windows with training data and inflate MCC by ~0.10 and AUC by ~0.04. Because 90% of UniProt "
    "palmitoylation annotations are inferred, and because inferred annotations are copied between "
    "orthologues, this is not a hypothetical risk: it is the default outcome of the splits used in most "
    "published work. We therefore release explicit cluster-held-out and species-held-out folds so that "
    "future methods can be compared on genuinely disjoint data."
))

rep("On those folds PalmDL reaches AUC 0.956 overall", (
    "On those folds PalmDL reaches AUC 0.944 overall and 0.914 on human sites. The best previously "
    "reported human AUC is 0.900 (GPS-Palm [8]); the best overall figure in the literature is 0.931 from "
    "a preprint using a site-level split (Deep-Palm [14]). Our human estimate is marginally above 0.900 "
    "(bootstrap P = 0.025) but its interval touches the published value, obtained under a protocol with "
    "zero shared windows. We therefore frame PalmDL as matching the published state of the art under a "
    "stricter protocol, with a small and uncertain advantage. This is not a head-to-head rerun: "
    "GPS-Palm's web server has been retired (HTTP 404 on 18 October 2026) and its standalone "
    "distribution is a GUI-only academic-licence package, so we could not score our test proteins with "
    "it. A true head-to-head on our released folds remains the single most important piece of missing "
    "evidence, and we provide the fold files to make it straightforward."
))

rep("The cross-species experiments split into two statements", (
    "The cross-species experiments split into two statements that are easy to conflate. When orthologues "
    "of the test species are available — the ordinary transfer setting, measured by the 6 × 6 matrix — "
    "mammalian models reach 0.85–0.95 on other mammals, effectively as good as within-species "
    "performance. When human-homologous clusters are removed, performance on human is 0.885. Most of the "
    "apparent generality of palmitoylation predictors therefore comes from orthologue coverage rather "
    "than from an abstract motif grammar. Bacterial lipidation does not belong in this comparison: after "
    "the label correction the bacterial positives are N-palmitoyl cysteines, a different modification "
    "with a different recognition logic, and we removed them rather than pretend that a transfer matrix "
    "over mixed chemistries measures the S-palmitoylation code."
))

rep("Two design choices follow from the error analysis.", (
    "Two design choices follow from the error analysis. First, the gap between inferred (AUC 0.947) and "
    "experimental sites (0.907) is an annotation-memory effect that any user should know about: the "
    "experimentally characterised sites are the ones the propagation machinery does not already cover. "
    "Second, reporting a fixed 0.5 cut-off is less useful than reporting FDR-controlled operating "
    "points, because the deployable question is \"how many predictions can I trust\", not \"what is my "
    "AUC\". At a 5% estimated FDR the model retains 34.5% of known sites at 90.3% precision, which is a "
    "realistic basis for a shortlist."
))

rep("Structure is a modest, real improvement rather than a cure.", (
    "Structure is a modest, real improvement rather than a cure. Three AlphaFold-derived numbers — "
    "pLDDT at the site, local pLDDT and C-alpha packing — improve the forest view from 0.890 to 0.911, "
    "consistently in four of five folds. They also point in the expected direction: palmitoylated "
    "cysteines are less ordered and more exposed than the cysteines that surround them. But the domain "
    "shift documented in §2.10 is not simply a burial-versus-exposure problem, and structure alone will "
    "not close it; a larger structured training set, enzyme conditioning or orthogonal experimental data "
    "will be needed."
))

# Limitations: replace the whole paragraph (it is long but starts uniquely)
rep("**Limitations.** (i) 93% of positives", (
    "**Limitations.** (i) 90% of positives are inferred rather than experimental; we stratify by evidence "
    "but cannot eliminate label noise, and the 449 experimental sites are too few for a well-powered "
    "experimental-only benchmark. (ii) Background negatives are proteins without annotation, not proteins "
    "proven not to be palmitoylated; some are certainly false negatives, which biases absolute metrics "
    "downwards while leaving model ranking intact. (iii) Our cluster definition uses 4-mer cosine "
    "similarity, an operational proxy: cosine 0.28–0.34 corresponded to median 71% global identity in a "
    "sample of 30 pairs, so remote homologues below ~50% identity may still cross folds; the species "
    "holdout is the more conservative experiment and is reported for that reason. (iv) We model only "
    "cysteine-centred windows and do not use predicted structure, DHHC-enzyme identity, or tissue "
    "expression. (v) The ensemble is calibrated at natural prevalence; performance at other prevalences "
    "requires recalibration. (vi) Independent, site-level public benchmarks do not exist: SwissPalm is "
    "behind an IP-verification gate, dbPTM and CysModDB no longer resolve, the MDD-Palm supplements "
    "contain no site lists, and GPS-Palm's server is retired. Both alternatives we could obtain fail a "
    "motif sanity check — the public RFCM-PALM mouse MS data (PWM AUC 0.53; PalmDL 0.52) and our own "
    "literature-mined set (PWM AUC 0.54; PalmDL 0.59; §2.10) — so neither is used as a benchmark. Their "
    "failure is itself informative: site-level assignment in cysteine-rich or transmembrane contexts is "
    "not yet reliable at scale, which is why the field's headline numbers remain independently "
    "unverified. (vii) The proteome-wide scan in this revision uses the forest view, because the ESM "
    "view could not be run at proteome scale on the available CPU hardware within a reasonable time; the "
    "released model and the prediction service still use the full ensemble, and the forest view is the "
    "stronger single view on hard negatives."
))

# ---- Methods ----
rep("Lipidation features were parsed with a regular expression", (
    "Lipidation features were parsed with a regular expression. A site was labelled positive only when "
    "the note was **S-palmitoyl cysteine** and the canonical residue at the annotated position was a "
    "cysteine; the first version instead matched any note containing \"palmitoyl\", which admitted "
    "2,381 N-/O-palmitoyl sites. Evidence codes were mapped to experimental (ECO:0000269, ECO:0007744, "
    "ECO:0007005, ECO:0006056) or inferred. Each site was represented by a 31-residue window centred on "
    "the cysteine and zero-padded with `X` at protein termini. Hard negatives were unannotated cysteines "
    "of positive-carrying proteins at least 16 residues from any palmitoylation site of any chemistry, "
    "capped at four per protein; background negatives were up to three cysteines from each human protein "
    "without any lipidation feature. Topology was derived from UniProt topological domains."
))

rep("Proteins carrying at least one positive were clustered", (
    "Every protein in the dataset — positive-carrying and background alike — was clustered by cosine "
    "similarity of its 4-mer count vector at a threshold of 0.30, and clusters sharing an identical "
    "31-residue window were merged, so no exact peptide occurs in two folds. Clusters were assigned to "
    "five folds with greedy balancing of positive counts (885 positives per fold). The species holdout "
    "used all human positives and human same-protein negatives as test; background negatives are human by "
    "construction and were excluded from both sides. Training used non-human sites after removing (i) "
    "every cluster containing a human protein and (ii) every annotation whose evidence referenced a human "
    "accession (`ECO:0000250|UniProtKB:<human>`). This leaves 1,480 training positives and 557 test "
    "positives."
))

rep("For each outer fold, 15% of the training clusters", (
    "For each outer fold, 15% of the training clusters were held out as an inner validation set for early "
    "stopping and threshold selection. Background negatives were subsampled to 1× the number of positives "
    "in training; test folds retained their natural composition. Thresholds were chosen to maximise MCC on "
    "the inner validation set, and a second threshold was chosen on a class-balanced subsample of the "
    "validation set for balanced metrics. Ensemble members were rank-transformed within fold and averaged "
    "with equal weights; a Platt scaling model was fitted on the inner validation scores and applied to "
    "test scores. Member subsets were compared by pooled validation AUC, never by test performance; the "
    "selected two-view combination achieved 0.9487, against 0.9428 for the forest, 0.9328 for ESM-2 and "
    "0.9420 for the three-view variant."
))

rep("Eight species with at least 85 annotated sites", (
    "Six species with at least 85 S-palmitoyl cysteine sites were used: human, mouse, rat, cattle, pig and "
    "*Arabidopsis thaliana*. For each source species a PalmDL model (random forest plus ESM head) was "
    "trained on that species' positives, its same-protein negatives and an equal number of human "
    "background negatives. For the diagonal, 20% of the source species' protein clusters were held out "
    "and the model retrained on the remaining 80%; off-diagonal entries were predicted by the model "
    "trained on all source rows, after removing target rows whose cluster appears in training. The ESM "
    "head used 15% of the remaining clusters for early stopping (30 epochs, patience 6). AUC was computed "
    "on each target species' positives and same-protein negatives."
))

rep("For each target FDR f, the validation folds were sorted", (
    "For each target FDR f and each fold separately, the validation rows were sorted by predicted "
    "probability and the largest threshold whose cumulative false-discovery fraction (FP / (TP + FP)) "
    "was ≤ f was selected. That fold's threshold was applied to that fold's outer test rows, and "
    "precision, recall and the number of predictions were pooled across folds. Thresholds are therefore "
    "never informed by other folds or by the test data."
))

rep("Three hundred randomly chosen human positive windows", (
    "Three hundred randomly chosen human positive windows were scored with the forest view. Each of the "
    "31 positions was mutated to each of the 20 amino acids in turn (620 predictions per site), and the "
    "mean change in predicted probability relative to the wild-type window was recorded. The same "
    "procedure exposed the label-filter bug described in §2.9, because mutating the central cysteine to "
    "lysine increased the score under the contaminated labels."
))

rep("All 262,286 cysteines of the 20,431 reviewed human sequences", (
    "All 262,286 cysteines of the reviewed human sequences in the local UniProt tables were scored in "
    "chunks of 5,000 with the forest view of the released model, calibrated on the same out-of-fold "
    "validation rows as the ensemble. The ESM view was not run at this scale in the revision because the "
    "CPU-only language model could not be evaluated over 262,286 windows within a reasonable time. For "
    "the FDR estimate we drew 20,000 random target windows and shuffled their 30 flanking residues while "
    "keeping the central cysteine fixed; this destroys the motif but preserves the candidate class, "
    "whereas mutating the centre (C→A) leaves a motif-shaped window and yields a badly conservative null "
    "because the model ignores the identity of the central residue. FDR at a threshold was estimated as "
    "(number of decoys at or above the threshold) × (targets / decoys) divided by (number of targets at "
    "or above the threshold), and q-values were obtained by taking the running minimum over lower "
    "thresholds. Because the decoy tail is thin, the resulting FDR should be read as optimistic, and the "
    "tiers are reported as candidate rankings rather than as proven error rates."
))

rep("Human proteins with at least one annotated site (316)", (
    "Human proteins with at least one annotated site and a random sample of background human proteins "
    "were queried through the AlphaFold DB prediction API, which returns the current model URL; models "
    "were retrieved for chains shorter than 700 residues, covering 1,950 human sites (463 positives in "
    "476 proteins). Per site we recorded the C-alpha pLDDT, the mean pLDDT of residues within ±5 "
    "positions, and the number of C-alpha atoms within 10 Å of the site C-alpha excluding sequence "
    "neighbours closer than four positions. The three features were standardised and appended to the base "
    "feature matrix; the same random forest and the same cluster-disjoint folds were used with and "
    "without them, so the comparison is paired at the fold level. The prediction service released with "
    "the paper exposes the sequence path by default."
))

# ---- Figure legends ----
rep("**Figure 1. Dataset composition.**", (
    "**Figure 1. Dataset composition.** (a) Organisms contributing the largest numbers of annotated "
    "S-palmitoyl cysteine sites (human and mouse highlighted; remaining 784 organisms pooled). (b) "
    "Evidence tier of positive sites (log scale): 449 experimentally supported versus 3,975 inferred, "
    "mostly by similarity. (c) Annotated topology of positive sites."
))
rep("**Figure 2. Split protocol determines measured performance.**", (
    "**Figure 2. Split protocol determines measured performance.** A single random forest was trained and "
    "evaluated under four protocols. (a) AUC and class-balanced MCC. (b) Percentage of test positives "
    "whose exact 31-residue window is present among training positives. Random site-level and "
    "protein-level splits leak 41–45% of windows; cluster and species splits leak none, for positives and "
    "negatives alike."
))
rep("**Figure 3. Discrimination under cluster-held-out cross-validation.**", (
    "**Figure 3. Discrimination under cluster-held-out cross-validation.** (a) ROC and (b) "
    "precision–recall curves for the PWM baseline, logistic regression, random forest, CNN, ESM-2 view and "
    "the PalmDL ensemble (4,424 positives, 67,266 sites)."
))
rep("**Figure 7. The palmitoylation code is eukaryotic.**", (
    "**Figure 7. The palmitoylation code is shared across mammals and partly with plants.** AUC of "
    "single-species PalmDL models, rows = training species, columns = test species, for the six species "
    "with at least 85 S-palmitoyl cysteine sites. Mammalian models transfer almost losslessly to other "
    "mammals; the plant model transfers only partially to mammals. Diagonal entries come from a "
    "protein-cluster-disjoint 20% holdout of the same species."
))
rep("**Figure 9. Proteome-wide scan and FDR-controlled shortlist.**", (
    "**Figure 9. Proteome-wide scan and candidate shortlist.** (a) Score distributions of all 262,286 "
    "human cysteines and of 20,000 shuffled-flank decoys (log scale). (b) Estimated FDR as a function of "
    "shortlist size, with the five reported tiers marked. (c) Composition of each tier: known UniProt "
    "sites versus novel cysteines."
))
rep("**Table 1.** Cluster-held-out cross-validation", (
    "**Table 1.** Cluster-held-out cross-validation (all species, 5 folds, 67,266 sites). AUC, AUPR, MCC "
    "and balanced MCC with 95% bootstrap CIs for PalmDL."
))
rep("**Table 5.** Cross-species transfer matrix", (
    "**Table 5.** Cross-species transfer matrix (AUC) for six species."
))
rep("**Table 8.** Proteome-wide scan", (
    "**Table 8.** Proteome-wide scan of the reviewed human proteome at five candidate tiers, with "
    "recovery of known sites and the fraction of novel predictions."
))

MD.write_text("\n".join(lines), encoding="utf-8")
print("done; remaining placeholders:", sum(1 for ln in lines if "{{" in ln))
