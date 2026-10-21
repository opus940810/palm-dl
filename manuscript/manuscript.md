# PalmDL: a leakage-aware deep learning framework and benchmark for protein S-palmitoylation site prediction

**Running title:** Leakage-aware deep learning for S-palmitoylation sites

**Authors:** [Author 1]¹, [Author 2]¹, [Corresponding author]¹\*

¹ [Affiliation], [City], [Country]

Correspondence: [email]

---

## Abstract

**Motivation.** S-palmitoylation is a reversible lipid modification that controls protein trafficking, signalling and membrane organisation, and its dysregulation underlies cancer, neurodegeneration and immune disease. Experimental mapping remains low-throughput, so dozens of computational predictors have been published. Yet the field has no common, homology-controlled benchmark: most tools are trained and tested on random or protein-level splits in which near-identical peptides, propagated annotations and orthologues cross the train/test boundary, and no predictor has been evaluated on a prospectively held-out species.

**Results.** We assembled a versioned, evidence-stratified dataset from UniProt release 2026_03 comprising 67,266 labelled cysteine sites (4,424 S-palmitoyl cysteine sites, of which 455 are experimentally supported; 786 organisms; 557 human sites in 312 proteins) with explicitly separated hard (same-protein) and background negatives. Restricting positives to S-palmitoyl cysteine excluded 2,381 N-/O-palmitoyl annotations that a permissive text filter had admitted, and the corrected folds share no identical 31-residue window between training and test. Using a single model across four protocols we show that conventional benchmarks leak severely: 41.1% of test positives share an identical window with a training positive under a random site-level split and 45.2% under a protein-level split, inflating AUC from 0.944 to 0.985 and MCC from 0.720 to 0.810. PalmDL, a two-view ensemble of ESM-2 embeddings and motif/physicochemical features whose member combination is chosen inside each outer fold, reaches AUC 0.938 (cluster-bootstrap 95% CI 0.917–0.954), AUPR 0.697 and MCC 0.590 under cluster-held-out cross-validation; on held-out human sites it reaches 0.910 (0.890–0.929) and 0.885 on a strict species holdout. It improves on each single view (ΔAUC +0.008 over ESM-2, +0.007 over the forest, whose interval includes zero). Probabilities are Platt-calibrated (mean fold ECE 0.024; pooled ECE 0.018). A six-species transfer matrix shows heterogeneous transfer rather than a universal code: most mammalian pairs transfer at 0.83–0.95, the plant model reaches 0.74–0.78 on mammals, and bacterial lipidation is excluded as a different modification. In-silico mutagenesis recovers the cysteine-cluster rule, and empirical validation-target operating points support candidate shortlists. A literature-mined external set of 136 sites from 91 papers remains motif-atypical (within-protein AUC 0.585), which we report as a domain-shift limitation rather than as validation.

**Availability.** All code, processed data, fold assignments and trained-model recipes are provided; the pipeline is reproducible end-to-end from public UniProt queries.

**Keywords:** S-palmitoylation; post-translational modification; deep learning; protein language model; benchmark leakage; cross-species generalisation

---

## 1. Introduction

S-palmitoylation is the reversible attachment of a 16-carbon fatty acid to cysteine thiols. Unlike most lipid modifications it is enzymatically reversible, and the competing actions of ~23 DHHC acyltransferases and several thioesterases make it a dynamic switch that partitions proteins between membrane compartments, controls synaptic transmission, and shapes oncogenic signalling [1–3]. The functional importance of the modification is matched by the difficulty of measuring it: acyl-biotin exchange, click-chemistry labelling and mass spectrometry can map sites, but each is biased by abundance, accessibility and labelling efficiency, and proteome-wide site-level coverage remains sparse.

Computational prediction has therefore been attractive for almost two decades. CSS-Palm 1.0 introduced clustering-based scoring with BLOSUM62 [4], CSS-Palm 2.0 added matrix mutation and a homology-reduced benchmark [5], and SVM-based and motif-based methods followed, including NBA-Palm [6] and MDD-Palm [7]. Deep learning arrived with GPS-Palm, a parallel CNN over sequence and structure-derived features that reported a human AUC of 0.900 [8], and pCysMod, a multi-PTM cysteine model reporting AUC 0.807 for palmitoylation [9]. More recent entries include RFCM-PALM for the mouse synaptic proteome [10], the SAPPTree structure-based random forest (preprint) and Deep-Palm, a stacking ensemble that reports AUC 0.931 on a CD-HIT 60% site-level 80/20 split (preprint [14]). CSS-Palm 3.0 and 4.0 were distributed as standalone software but have no peer-reviewed description, and their web servers have been retired.

Three problems limit the comparability and credibility of this literature. First, negatives are almost always "all other cysteines" of the same proteins or random cysteines from unannotated proteins; the ratio and difficulty of negatives differ between studies, and accuracy or AUC computed on easy background negatives is not comparable across papers. Second, split protocols are rarely protein- or family-disjoint. Palmitoylation annotations are heavily propagated by similarity: in UniProt, 90% of S-palmitoyl cysteine annotations are inferred rather than experimentally determined, and orthologues across species carry identical 31-residue windows. A random site-level split therefore places the same peptide in training and test, and even a protein-level split does so whenever annotations were propagated between species. Third, no palmitoylation predictor has been tested on a species that was entirely absent from training, so it is unknown whether learned motifs transfer or whether they memorise species- and family-specific sequence context.

Here we address all three. We build a versioned, evidence-stratified dataset from a single UniProt release and separate hard from background negatives. We quantify, with one fixed model, how much each splitting protocol inflates performance, and measure the fraction of test positives whose exact window appeared in training. We then develop PalmDL, a leakage-aware ensemble of a protein-language-model view (ESM-2) and a motif/physicochemical view (random forest), calibrated and evaluated under cluster-held-out cross-validation, human-only training, and a strict cross-species holdout. We report PalmDL's human AUC on a benchmark with zero shared windows alongside published values obtained on different data, and we release the split definitions so that future methods can be compared on identical, leakage-controlled folds.

---

## 2. Results

### 2.1 A versioned, evidence-stratified dataset

We queried UniProt release 2026_03 for all Swiss-Prot entries carrying at least one lipidation feature (9,089 entries) and kept only features whose note is **S-palmitoyl cysteine** on a cysteine residue. The first version of this filter matched any note containing "palmitoyl" and therefore admitted N-palmitoyl cysteine, N6-palmitoyl lysine and O-palmitoyl Ser/Thr; the corrected dataset contains **4,424 positive sites**, 455 supported by experimental evidence codes and 3,969 inferred (mostly by similarity, 89.7%), from 786 organisms. Human (557 sites in 312 proteins) and mouse (511) are the largest contributors (Figure 1a,b). Of the annotated positives, 2,445 are cytoplasmic and only 27 extracellular (Figure 1c).

Negatives were defined in two tiers. **Hard negatives** are unannotated cysteines of the same proteins whose ±15-residue windows do not overlap any palmitoylation annotation, whatever its chemistry (8,620 sites). **Background negatives** are cysteines sampled from 19,622 reviewed human proteins that carry no lipidation feature at all (54,222 sites). The final labelled set contains 67,266 cysteine sites from 21,568 proteins and 61,183 distinct 31-residue windows. Metrics remain high on the harder same-protein negatives alone (AUC 0.952; Table 3), so model ranking is not an artefact of easy background negatives.

### 2.2 Conventional splits leak, and the leakage is measurable

We trained the same random forest on the same features under four protocols and evaluated each on holdout sites that were never used for fitting or threshold selection (Table 2, Figure 2). Under a random site-level split, **41.1% of test positives had an identical 31-residue window among the training positives**, and AUC reached 0.980 (MCC 0.843). A protein-level split, the most common precaution in the field, still leaked **45.2%** of test windows and gave AUC 0.985 (MCC 0.810). Cluster-level splitting removed the overlap entirely — **0% of any test row, positives and background negatives alike, shares a window with training** — and reduced AUC to 0.944 and MCC to 0.720. Holding out an entire species gave AUC 0.885.

This audit has two consequences. First, headline numbers from site-level or protein-level splits cannot be compared with numbers from homology-controlled splits; a model reporting AUC > 0.97 is likely measuring peptide memorisation. Second, exact-window deduplication must cover the negative pool as well: in the first version of our splits, background proteins were not clustered and 751 windows (2,012 rows, almost all background negatives) crossed folds. Clustering every protein in the dataset and merging groups that share a window removes that residue. We therefore report every result in this paper under cluster-level splits, and additionally under a species holdout.

### 2.3 Model development: two complementary views

We compared five model families under identical cluster-held-out folds (Table 1, Figure 3). A CSS-Palm-style position-specific log-odds matrix reached AUC 0.799 and MCC 0.344. Logistic regression on the flattened one-hot window reached AUC 0.829. A multi-kernel CNN reached AUC 0.854 and MCC 0.439. A random forest on one-hot and physicochemical channels reached AUC 0.931 and MCC 0.595. An MLP on frozen ESM-2 (35M) embeddings reached AUC 0.931 and MCC 0.612. After the label correction the two strongest views are therefore tied on discrimination, with ESM-2 better at fixing a decision threshold.

The two views make different errors: the forest captures exact positional chemistry and local motif syntax, whereas ESM-2 embeddings capture evolutionary and structural context. Member selection is now **nested inside each outer fold**: the 31 non-empty subsets are ranked by that fold's validation AUC, the winner is calibrated on the same validation rows, and only then applied to that fold's test rows. Four folds select RF+ESM and one selects LR+RF+CNN+ESM. The resulting ensemble reaches **AUC 0.938 (cluster-bootstrap 95% CI 0.917–0.954), AUPR 0.697, MCC 0.590 and balanced MCC 0.747** (Table 1). Paired cluster bootstrap gives ΔAUC = +0.008 (0.002–0.014) over ESM-2 and +0.007 (−0.004–0.018) over the forest, so the advantage over the forest is not resolved by these folds. For comparison, the previous globally selected combination scores 0.944 (0.926–0.957 cluster), which shows that pooling validation folds to choose one combination is mildly optimistic; we report the nested estimate as primary.

### 2.4 Human evaluation and the cost of a strict split

Because most published predictors are evaluated on human sites, we repeated the analysis on human-only test folds (Table 3, Figure 4). With multi-species training and human test folds, nested PalmDL reaches **AUC 0.910 (0.890–0.929), AUPR 0.294, MCC 0.297 and balanced MCC 0.662**. These values reproduce the published human AUC of 0.900 only descriptively: no paired test is possible across different test sets, and the within-study bootstrap probability that our estimate falls below 0.900 is 0.08. We therefore report the number and leave the comparison to context. Training on human data alone is selected as RF+ESM in four folds and CNN+ESM in one (mean 0.901). On the strict species holdout the random forest reaches AUC 0.885 and balanced MCC 0.614. The 0.025 difference between the human CV estimate and the species holdout uses different model families and negative compositions, so it is descriptive rather than a clean decomposition of a species-specific cost.

### 2.5 Calibration, hard negatives and interpretability

Probabilities from the calibrated ensemble are reliable at natural prevalence (mean fold ECE 0.024, pooled 10-bin ECE 0.018; Figure 5), whereas the uncalibrated forest has ECE ≈ 0.12 and over-confident tails. On the hard-negative benchmark (positives plus same-protein cysteines only, 1:1.95 ratio), PalmDL reaches AUC 0.952 and MCC 0.664, compared with 0.811/0.414 for the PWM baseline (Table 3), so the advantage is not an artefact of easy background negatives.

The learned motif agrees with known biochemistry (Figure 6). Positive windows are enriched in hydrophobic and aromatic residues at positions −3 to −1 relative to the cysteine (leucine, phenylalanine, tryptophan, tyrosine and methionine), depleted of glycine, proline and acidic residues in the same interval, and enriched in basic residues at positions +3 to +4. This is consistent with the DHHC active-site geometry, in which the cysteine is presented in a hydrophobic groove, and with the known preference for cysteine clusters and juxtamembrane positions.

### 2.6 Position relative to published predictors

Table 4 places our results next to the numbers reported by previous methods. Three patterns are visible. First, classical MCC values (0.45–0.67) were obtained under protocols that either reuse the training literature or do not hold out proteins, and the only classical tool with an independent test (CSS-Palm 2.0) dropped from MCC 0.671 to 0.523 when it did so. Second, the deep-learning systems with published AUCs (GPS-Palm 0.900, pCysMod 0.807) report neither cross-species nor explicit protein-disjoint evaluation. Third, the highest published number (Deep-Palm 0.931) comes from a site-level split, precisely the protocol that our audit shows to leak 41–45% of windows. PalmDL's 0.938 is obtained with zero shared windows under nested member selection, and its human estimate, 0.910, is reported for context rather than as a paired superiority claim, because the published values come from different data, negatives and splits. Because we could not rerun GPS-Palm or the retired classical tools on our folds, Table 4 compares protocols as much as algorithms; the released fold files remove this obstacle for future work.

**Table 4. Comparison with published S-palmitoylation predictors.**

| Predictor | Year | Evaluation protocol | Reported performance | Protein/species-disjoint test |
|---|---|---|---|---|
| CSS-Palm 1.0 | 2006 | jackknife, 210 sites | MCC 0.588 | no |
| NBA-Palm | 2006 | re-evaluated on CSS-Palm data | MCC 0.610 | no |
| CSS-Palm 2.0 | 2008 | leave-one-out | MCC 0.671 | yes: MCC 0.523 on 53 later sites |
| MDD-Palm | 2017 | 5-fold, 710 sites | MCC 0.45 | partial (Sn/Sp 0.69) |
| GPS-Palm | 2021 | internal, 3,098 sites | human AUC 0.900 | not reported |
| pCysMod | 2021 | internal | AUC 0.807 | not reported |
| RFCM-PALM | 2021 | mouse synaptic holdout | MCC 0.532 | mouse only |
| Deep-Palm (preprint) | 2026 | site-level 80/20, CD-HIT 60% | AUC 0.931 | no |
| **PalmDL** | this work | cluster-held-out CV, 4,424 sites | **AUC 0.938, MCC 0.590** | yes (0% shared windows) |
| **PalmDL** | this work | human test folds | **AUC 0.910, balanced MCC 0.662** | yes (clusters) |
| **PalmDL** | this work | species holdout | AUC 0.885, balanced MCC 0.614 | yes (strictest) |

### 2.7 Cross-species transfer is heterogeneous

The species holdout in §2.4 removes every human-homologous cluster and therefore measures transfer to genuinely novel families. To separate that from ordinary cross-species transfer we trained one model per species on its positives and same-protein negatives, holding out clusters for the diagonal and for early stopping, and evaluated every model on every species with target rows restricted to positives and same-protein negatives (Figure 7, Table 5). Transfer is **heterogeneous**, not uniformly lossless. Most mammalian pairs fall between 0.83 and 0.95; the arithmetic mean of the 18 observable mammalian off-diagonal cells is 0.901. Pig→Human (0.834) is the clear exception, 0.084 below the human diagonal. Pig's own diagonal is empty (fewer than 20 eligible positives), so no diagonal comparison is made for that row. The plant model reaches 0.739–0.785 on mammals, whereas mammalian models reach 0.96–0.98 on *A. thaliana*. Several cells have too few eligible positives after cluster removal and are left empty rather than estimated from a handful of sites. The matrix therefore describes transfer among these sampled species and negative-sampling protocols; it does not isolate the causal contribution of orthologue coverage.

Two consequences follow. First, claims of a single shared eukaryotic code are not supported by this design: the pig exception and the plant asymmetry are as informative as the mammalian mean. For the human source only, the early-stopping split also contained 6,298 background negatives (the other sources use same-protein negatives), so the human row's model-selection condition is not fully matched; this is disclosed rather than corrected in this revision. Second, the strict human holdout (0.885) and the human CV estimate (0.910) differ by only 0.025 and use different model families and negative compositions, so they should not be read as a clean species-specific penalty. A paired design with matched negatives and a frozen target set would be needed to separate orthologue coverage from motif generalisation.


**Table 5. Cross-species transfer matrix (AUC).**

| train \ test | Human | Mouse | Rat | Cattle | A. thaliana | Pig |
|---|---|---|---|---|---|---|
| Human | 0.917 | 0.928 | 0.910 | 0.948 | 0.980 | — |
| Mouse | 0.899 | 0.877 | 0.880 | 0.883 | 0.984 | — |
| Rat | 0.890 | 0.914 | 0.928 | 0.945 | 0.984 | 0.940 |
| Cattle | 0.899 | 0.903 | 0.917 | 0.880 | 0.979 | 0.937 |
| *A. thaliana* | 0.739 | 0.766 | 0.781 | 0.775 | — | 0.785 |
| Pig | 0.834 | 0.854 | 0.858 | 0.880 | 0.964 | — |

One model per source species is trained on its positives and same-protein negatives (no human background); 20% of source clusters are held out for the diagonal and another 15% for early stopping. For the human source the early-stopping split also contained 6,298 background negatives, unlike the other rows; the 18-cell mammalian off-diagonal mean is 0.901. Target cells use only positives and same-protein negatives, and rows whose cluster appears in either the fit or the validation split are removed. Off-diagonal cells therefore share the model of their row, as stated in Methods. Dashes mark cells with fewer than 20 eligible test positives (Human→Pig 13, Mouse→Pig 19, and the *A. thaliana* and Pig diagonals with 15 and 16).

Diagonal entries are from a protein-cluster-disjoint 20% holdout of the same species; off-diagonal entries use a model trained on all source rows, with target rows in the training clusters removed. Pig and A. thaliana have few held-out positives, so their diagonal is correspondingly noisy.


### 2.8 Error structure and validation-target operating points

Because a predictor is only useful if its errors are predictable, we stratified the held-out predictions (Table 6). Discrimination remains better for inferred than for experimentally supported sites (pooled AUC 0.941 versus 0.902; each tier uses its positives plus all 62,842 negatives). Evidence tiers now follow one rule: a site counts as experimental when **any** evidence clause is experimental (ECO:0000269, 0007744, 0007005 or 0006056), which gives 455 experimental and 3,969 inferred sites; six mixed-clause sites were previously misclassified by a comma-only split and are now aligned with Table 6. Pooled performance is stable across annotated topologies (cytoplasmic 0.858, extracellular 0.828, topology unknown 0.900); the extracellular estimate draws on three of five folds with both classes and the lumenal class (14 positives) is not interpretable. These pooled values differ from the unweighted fold means in Table 1, which are the primary estimates. Topology is itself an input feature, so these strata describe where the model is accurate, not an independent validation of the topology signal.

For practical use we derived, **within each fold**, the lowest threshold whose complete score groups satisfy an empirical validation-target operating point, and evaluated it on that fold's outer test rows (Figure 8b). The nominal targets are not achieved out of sample, and we report both sides: at a 1% validation target the held-out precision is 0.953, recall 0.201 and observed false-discovery proportion 4.7%; at 5%, 0.900/0.331/10.0%; at 10%, 0.849/0.398/15.1%; at 20%, 0.759/0.490/24.1%. These are empirically calibrated operating points, not a demonstrated FDR-controlling procedure; the gap between validation target and held-out FDP is the honest measure of their transferability.

**Table 6. Error structure of the held-out ensemble predictions.**

| Stratum | n positives | AUC | AUPR |
|---|---|---|---|
| Experimentally supported sites | 455 | 0.902 | 0.176 |
| Inferred sites | 3,969 | 0.941 | 0.664 |
| Cytoplasmic | 2,445 | 0.858 | 0.808 |
| Extracellular | 27 | 0.828 | 0.361 |
| Topology unknown | 1,801 | 0.900 | 0.512 |

### 2.9 What the model uses: cysteine clusters, charge exclusion, and a label audit

In-silico mutagenesis of 300 human positive windows (every position mutated to every amino acid, scored with the forest view; Figure 8a) shows that substitutions at the central position are neutral, as they must be: every candidate is a cysteine, so the centre carries no discriminative information. The informative positions are the flanks. Introducing a cysteine at +1 increases the predicted probability by 0.089 and at −1 by 0.083, recovering the known preference for cysteine clusters; introducing glutamate or aspartate at −1 decreases it by 0.06–0.09, recovering the exclusion of acidic residues near the modification site.

The same audit exposed a bug in **our own label filter**, not in UniProt. The first version of the dataset matched any lipidation note containing "palmitoyl", which admitted N-palmitoyl cysteine (2,315 sites), N6-palmitoyl lysine (34), N-palmitoyl glycine (18) and O-palmitoyl serine/threonine (14) — 35% of the original positives. These are legitimate annotations of different chemical modifications; they were mis-classified by us, and the earlier claim that they reflected UniProt annotation errors was wrong and is withdrawn. The corrected dataset requires an S-palmitoyl cysteine note on a cysteine residue and contains no non-cysteine positives. The label correction changed every downstream number in this paper; the audit trail (per-site notes before and after) is released with the data, and we report the corrected results throughout.

### 2.10 Literature mining cannot yet supply an independent benchmark

Finally we asked whether the missing external benchmark can be manufactured from the literature. We searched Europe PMC for open-access articles on palmitoylation published between 2020 and 2026 (9,062 hits), fetched 992 full texts, and extracted every sentence that mentions palmitoylation together with a cysteine position. Protein tokens were resolved against 28,711 locally cached reviewed sequences, and a mention was kept only when the residue at the stated position is a cysteine: **136 sites in 109 proteins from 91 papers** (Table 7). Only 2.2% of these sites are annotated in UniProt 2026_03, so the set lies almost entirely outside the training labels and is the closest thing to an independent site-level benchmark currently obtainable.

It is not a clean benchmark. The motif that the model learned is barely enriched in this set (PWM AUC 0.537 on the 133 scorable sites), and PalmDL separates the literature sites from same-protein cysteines only weakly (weighted within-protein AUC 0.585; pooled novel-site AUC 0.578). Manual inspection explains why: genuine sites are mixed with extraction errors, and many genuine sites are transmembrane-proximal or otherwise motif-atypical cysteines that an annotation-trained model has least often seen. A low PWM score does not by itself prove that a literature site is wrong; it shows that this class lies outside the motif distribution the model was trained on.

To make sure this was not an extraction artefact, we curated all 136 sites by reading every sentence. Ninety-nine mentions could be assigned unambiguously to a protein and a cysteine — either because the sentence states the site explicitly or because a cysteine-to-serine/alanine mutant lost palmitoylation — while 37 were rejected as sites of a different protein in the same sentence, catalytic-site mutants of DHHC enzymes, predicted-only sites, negative results, or table rows with mismatched proteins. The curated subset also has little overlap with the candidate tiers (0/1/2 at 5%/20%/30%), consistent with the full set; extraction noise alone is not a sufficient explanation, but the curated set has no independent AUC of its own.

The curated sites are almost absent from the candidate list: **none of the 99 appears at the 5% tier, one at the 20% tier and two at the 30% tier** (the same counts as the full 136-site set). For PalmDL, the honest conclusion is a domain-shift limitation: it reproduces and extends the UniProt consensus, but does not recognise the motif-atypical, transmembrane-proximal sites that mechanistic papers tend to characterise.

We therefore report this as an audit of the resource landscape rather than as a validation. It quantifies the statement in the Introduction for our own model: its headline numbers have not been checked on data that are simultaneously independent, site-level and motif-consistent, and automated literature mining does not yet close that gap. A low motif score does not establish that a literature site is wrong; it shows that the site lies outside the motif distribution the model learned.

**Table 7. Automated literature mining as an external benchmark.**

| Filter | sites | PWM AUC | PalmDL AUC |
|---|---|---|---|
| All verified sites | 133 | 0.537 | 0.585† |

†Weighted within-protein AUC over 104 protein-level tests; the pooled novel-site AUC against same-protein cysteines is 0.578 and the filter-table value for the same 133 sites is 0.588. Comparisons are against other cysteines of the same proteins. The manually curated subset (99 unambiguous mentions in 82 proteins) is used for the overlap analysis, not as a separate filter row.

### 2.11 A proteome-wide candidate list with empirical q-values

As a forward-looking counterpart to the retrospective benchmarks, we applied the forest view of the released model to every cysteine of the reviewed human proteome — **262,286 sites in 19,814 proteins** — and estimated empirical q-values against 20,000 shuffled-flank decoys that keep the central cysteine but destroy the motif. (A C→A decoy is inappropriate: §2.9 shows that the model is deliberately insensitive to the identity of the centre.) At the 1% tier the scan yields 215 sites in 134 proteins; at 5%, 301 sites in 199 proteins (248 already annotated, 53 novel); at 10%, 396; at 20%, 615; and at 30%, 1,010 (Table 8, Figure 9).

The 5% tier is biologically coherent in composition. 61.8% of its sites carry annotated cytoplasmic topology versus 5.6% of all scanned cysteines (11.1-fold). The highest-confidence examples at the 5% tier include endothelin receptor A C386 (q = 0.044); its neighbour C388 is in the 20% tier (q = 0.174), and the other candidates mentioned in earlier drafts (claudin-2 C185, mucolipin-3 C550) have q = 0.70 and 0.41 and are not top-tier examples. Topology is an input feature, so this enrichment describes the composition of the list rather than validating the predictions.

Four caveats bound the resource. First, the proteome-wide prior is small (557 annotated sites among 262,286 cysteines, ≈0.2%), so any strict operating point is intrinsically severe. Second, 82% of scanned cysteines have no topology annotation and are strongly depleted from the shortlist (38.2%), so the list is biased toward proteins with curated topology. Third, the decoy tail is thin — few shuffled decoys reach the top of the score distribution — so the FDR estimate is optimistic and the tiers are candidate rankings rather than proven error rates. Fourth, known-site recovery in the scan (44.5% at the 5% tier) is an in-sample upper bound because the released model was trained on all UniProt annotations; the corresponding out-of-fold recall is 33.1% (§2.8).

**Table 8. Proteome-wide scan of the reviewed human proteome (262,286 cysteines) at five empirical candidate tiers.**

| empirical tier | sites | proteins | known recovered | recall of known | novel | score threshold |
|---|---|---|---|---|---|---|
| 1% | 215 | 134 | 197 | 35.4% | 18 | 0.979 |
| 5% | 301 | 199 | 248 | 44.5% | 53 | 0.959 |
| 10% | 396 | 269 | 276 | 49.6% | 120 | 0.938 |
| 20% | 615 | 416 | 334 | 60.0% | 281 | 0.894 |
| 30% | 1,010 | 701 | 380 | 68.2% | 630 | 0.830 |

"Known" refers to the 557 S-palmitoyl cysteine sites annotated in UniProt 2026_03 for human proteins; "novel" sites are unannotated cysteines. FDR is estimated by target–decoy comparison against 20,000 shuffled-flank decoys scored with the forest view; the decoy tail is thin, so these tiers are candidate rankings rather than proven error rates.

### 2.12 A structure-aware prototype: useful inside the annotation distribution, not outside it

Sequence-only models cannot see whether a cysteine is buried in a folded domain or exposed on a membrane-proximal surface, which matters because the acyl chain must be accommodated. We fetched AlphaFold models for the human proteins carrying annotated sites and computed three features per cysteine: pLDDT at the site, mean pLDDT within ±5 residues, and C-alpha packing density (number of C-alpha atoms within 10 Å, excluding sequence neighbours). 1,950 human sites (463 positives in 476 proteins) were covered.

Adding the features to the random-forest view lifts AUC from 0.890 to 0.911 on the structured subset (Δ +0.021, positive in 4 of 5 folds; Table 9). The effect is modest because the language-model view already captures much of the structural context, and because only annotated human sites have models. We therefore report structure as a paired analysis rather than as part of the released sequence model.

The direction of the structural effect is biologically coherent: palmitoylated cysteines have lower pLDDT (66.2 versus 81.1) and lower packing (4.9 versus 10.8 C-alpha neighbours), i.e. they sit in less ordered, more exposed segments where the acyl chain can be accommodated and where DHHC enzymes act.

No structural model was run for the literature sites in this revision, so whether structure would help there is untested. The earlier claim that structure cannot close the domain shift is therefore withdrawn: the paired comparison only establishes a modest gain where annotated human sites have AlphaFold models.

**Table 9. Structure-aware variants on human sites (paired cluster-disjoint folds).**

| Variant | Comparison | AUC | Δ |
|---|---|---|---|
| Forest + structure | 0.890 → 0.911 (1,950 sites) | 0.911 | +0.021 (4/5 folds) |

| Feature | Positive sites (mean) | Negative sites (mean) |
|---|---|---|
| pLDDT at cysteine | 66.2 | 81.1 |
| C-alpha neighbours within 10 Å | 4.9 | 10.8 |

---

## 3. Discussion

We set out to build a palmitoylation predictor that could credibly claim state-of-the-art performance. Doing so required first showing that the field's benchmarks overstate performance. With a single fixed model we measured that random and protein-level splits share 41–45% of test-positive windows with training data and inflate MCC by ~0.10 and AUC by ~0.04. Because 90% of UniProt palmitoylation annotations are inferred, and because inferred annotations are copied between orthologues, this is not a hypothetical risk: it is the default outcome of the splits used in most published work. We therefore release explicit cluster-held-out and species-held-out folds so that future methods can be compared on genuinely disjoint data.

On those folds nested PalmDL reaches AUC 0.938 overall and 0.910 on human sites. The best previously reported human AUC is 0.900 (GPS-Palm [8]); the best overall figure in the literature is 0.931 from a preprint using a site-level split (Deep-Palm [14]). These are different test sets and negative schemes, so the numbers are context rather than a paired comparison; the within-study bootstrap probability that our human estimate falls below 0.900 is 0.08. This is not a head-to-head rerun: GPS-Palm's web server has been retired (HTTP 404 on 18 October 2026) and its standalone distribution is a GUI-only academic-licence package. A true head-to-head on our released folds remains the single most important piece of missing evidence, and we provide the fold files to make it straightforward.

The cross-species experiments are deliberately descriptive. With orthologues available, most mammalian pairs transfer at 0.83–0.95 (18-cell mean 0.901), but Pig→Human is 0.834 and the plant→mammal direction is 0.739–0.785; without human-homologous clusters, human performance is 0.885. The matrix therefore shows heterogeneous transfer among sampled species and negative-sampling protocols, not a uniform eukaryotic code, and it does not isolate orthologue coverage as the cause. Bacterial lipidation is excluded rather than compared: after the label correction the bacterial positives are N-palmitoyl cysteines, a different modification.

Two design choices follow from the error analysis. First, the gap between inferred (AUC 0.947) and experimental sites (0.907) is an annotation-memory effect that any user should know about. Second, reporting a fixed 0.5 cut-off is less useful than reporting empirical operating points, provided their transfer is shown: at a 5% validation target the held-out precision is 0.900 and the observed false-discovery proportion 10.0%, so the nominal target is not a guarantee.

The proteome-wide scan turns that principle into a resource. At the 5% tier it returns 301 human cysteines, 53 of them unannotated, with an 11-fold enrichment for cytoplasmic topology and a top tier dominated by juxtamembrane cysteines of single-pass membrane proteins. The small absolute size is not a weakness of the model but a property of the problem: with only ~0.2% of human cysteines annotated, any procedure that controls FDR at the proteome scale must be conservative. We therefore report five tiers so that users can trade confirmation for discovery; the 20–30% tiers are the useful range for follow-up, and the decoy tail is thin enough that the tiers should be read as candidate rankings rather than proven error rates.

Structure is a modest, real improvement rather than a cure. Three AlphaFold-derived numbers — pLDDT at the site, local pLDDT and C-alpha packing — improve the forest view from 0.890 to 0.911, consistently in four of five folds. pLDDT is a model-confidence measure that correlates with disorder, not a direct experimental measure of flexibility. No structural model was available for the literature sites in this revision, so whether structure would help there remains untested rather than refuted.

**Limitations.** (i) 90% of positives are inferred rather than experimental; we stratify by evidence but cannot eliminate label noise, and the 455 experimental sites are too few for a well-powered experimental-only benchmark. (ii) Background negatives are proteins without annotation, not proteins proven not to be palmitoylated. (iii) The cluster definition uses 4-mer cosine similarity, an operational proxy: cosine 0.28–0.34 corresponded to median 71% global identity in 30 sampled pairs, so remote homologues below ~50% identity may still cross folds. (iv) We model only cysteine-centred windows and do not use predicted structure, DHHC-enzyme identity, or tissue expression. (v) The ensemble is calibrated at natural prevalence; performance at other prevalences requires recalibration. (vi) Independent, site-level public benchmarks do not exist: SwissPalm is behind an IP-verification gate, dbPTM and CysModDB no longer resolve, the MDD-Palm supplements contain no site lists, and GPS-Palm's server is retired. Both alternatives we could obtain fail a motif sanity check, so neither is used as a benchmark. (vii) The proteome-wide scan uses the forest view and a shuffled-flank null with a thin tail; the tiers are candidate rankings, not demonstrated error rates. (viii) The cross-species matrix samples heterogeneous negative sets and leaves four cells empty for lack of eligible positives; it does not support a causal orthologue-coverage claim. (ix) The deployed service averages view logits rather than within-fold ranks, which on the same folds gives AUC 0.940 versus 0.944 for the rank average under the global combination.

**Future work.** Structure-aware features from AlphaFold, enzyme-conditioned models that predict which DHHC isoform a site prefers, manually curated literature sets assembled with the extraction pipeline of §2.10, and a proteome-wide scan with target–decoy FDR estimation are the natural next steps. The released `score()` function makes a batch prediction interface straightforward to deploy.

---

## 4. Methods

### 4.1 Data acquisition

All sequences and annotations were retrieved from the UniProtKB REST API (release 2026_03) on 18 October 2026 using two queries: `(ft_lipid:*) AND (reviewed:true)` and `(reviewed:true) AND (organism_id:9606) AND (NOT ft_lipid:*)`. Fields retrieved: accession, sequence, organism, protein and gene names, lipidation features, signal peptide, topological domains and length. Queries and row counts are stored in `data/raw/download_meta.json`.

### 4.2 Site and negative definition

Lipidation features were parsed with a regular expression. A site was labelled positive only when the note was **S-palmitoyl cysteine** and the canonical residue was a cysteine. A site was classified as experimental when at least one ECO clause in its evidence string matched ECO:0000269, ECO:0007744, ECO:0007005 or ECO:0006056; clauses are separated by ',' or ';' in UniProt, and both separators are handled. This gives 455 experimental and 3,969 inferred sites, the counts used in Table 6. Each site was represented by a 31-residue window centred on the cysteine and zero-padded with `X`. Hard negatives were unannotated cysteines of positive-carrying proteins at least 16 residues from any palmitoylation site, capped at four per protein; background negatives were up to three cysteines from each human protein without any lipidation feature. Topology was derived from UniProt topological domains.

### 4.3 Leakage-controlled splits

Every protein in the dataset — positive-carrying and background alike — was clustered by cosine similarity of its 4-mer count vector at 0.30, and clusters sharing an identical 31-residue window were merged, so no exact peptide occurs in two folds. Clusters were assigned to five folds with greedy balancing of positive counts (885 positives per fold). The species holdout used human positives and human same-protein negatives as test; background negatives are human by construction and were excluded from both sides. Training used non-human sites after removing every cluster containing a human protein and every annotation propagated from a human accession, leaving 1,480 training and 557 test positives.

### 4.4 Features

Three feature blocks were computed per site. (i) A 23 × 31 one-hot/physicochemical tensor: 20 amino-acid channels plus Kyte–Doolittle hydropathy, normalised van der Waals volume and Vihinen flexibility. (ii) Eleven protein-context features: six topology indicators, an annotation-available flag, log protein length, relative position, relative distance to the nearest terminus, and a terminus-proximal flag. (iii) ESM-2 (esm2_t12_35M_UR50D) representations, computed as the mean over the 31 window tokens and the representation of the central cysteine (2 × 480 dimensions). ESM embeddings were extracted once and cached; no ESM fine-tuning was performed.

### 4.5 Models and training

The PWM baseline computes position-specific log-odds from training positives and negatives with a pseudocount of 0.5 and converts the summed score to a probability with a logistic link. Logistic regression and random forests (400 trees, minimum leaf size 3, balanced subsampling) used the flattened one-hot window plus context features. The CNN has four parallel 1-D convolutions (kernels 3, 5, 7, 11; 32 filters each), ReLU, global max pooling, batch normalisation, a 64-unit hidden layer with dropout 0.3, and a sigmoid output; it is trained with AdamW (learning rate 1 × 10⁻³, weight decay 1 × 10⁻⁴), class-weighted BCE, gradient clipping at 5 and early stopping on validation AUC (patience 8, maximum 40 epochs, batch size 128). The ESM head is a two-layer MLP (256 and 128 units, dropout 0.3) trained with the same schedule and a learning rate of 5 × 10⁻⁴.

For each outer fold, 15% of the training clusters were held out as inner validation. Background negatives were subsampled to 1× the number of positives in training; test folds retained their natural composition. Thresholds were chosen to maximise MCC on the inner validation set. Ensemble member selection, calibration and thresholding were performed **within the same outer fold**: the 31 non-empty subsets were ranked by that fold's validation AUC, the winner was Platt-calibrated on those rows, and the calibrated scores were applied to that fold's test rows. The globally selected combination is reported only as a sensitivity analysis. Discrimination is the unweighted mean of the five fold metrics; MCC and accuracy pool confusion counts at the per-fold thresholds. Two bootstrap schemes are reported: resampling sites and resampling protein clusters within folds.

### 4.6 Evaluation

Discrimination was measured by AUC and average precision; threshold-dependent performance by MCC, accuracy, sensitivity and specificity at the validation-derived threshold, reported both at natural prevalence and on class-balanced subsamples. Calibration was measured by expected calibration error with ten bins. Confidence intervals were obtained by bootstrapping sites (2,000 replicates) and, because sites within a protein cluster are not independent, by resampling protein clusters within folds (1,000 replicates); the wider cluster interval is the one quoted in the abstract. Model comparisons used paired cluster bootstrap on identical rows. All experiments used Python 3.12, PyTorch 2.9 (CPU), scikit-learn 1.6 and Biopython 1.88.

### 4.7 Cross-species transfer matrix

Six species with at least 85 S-palmitoyl cysteine sites were used: human, mouse, rat, cattle, pig and *Arabidopsis thaliana*. For each source species one model (random forest plus ESM head) was trained on that species' positives and its same-protein negatives; human background negatives were excluded. Twenty percent of the source clusters were held out for the diagonal and another 15% for early stopping. Off-diagonal cells use the same model as their row. Target rows were restricted to positives and same-protein negatives and removed if their cluster appeared in either the fit or the validation split. Cells with fewer than 20 eligible test positives are reported as empty rather than estimated.

### 4.8 Empirical validation-target operating points

For each empirical target f and each fold separately, validation rows were sorted by predicted probability and the lowest threshold was chosen whose **complete score groups** satisfied the cumulative false-discovery constraint (FP / (TP + FP) ≤ f); splitting a tied group is not allowed. The threshold was applied to that fold's outer test rows and precision, recall, predicted counts and the observed test false-discovery proportion were pooled across folds. Validation targets are nominal: they are not a guarantee of out-of-sample FDR control, and both values are reported.

### 4.9 In-silico mutagenesis and label audit

Three hundred randomly chosen human positive windows were scored with the forest view. Each of the 31 positions was mutated to each of the 20 amino acids in turn (620 predictions per site), and the mean change in predicted probability relative to the wild-type window was recorded. The same procedure exposed the label-filter bug described in §2.9, because mutating the central cysteine to lysine increased the score under the contaminated labels.

### 4.10 Literature mining

Europe PMC was queried for `(palmitoylation OR "S-palmitoylation" OR palmitoylated) AND OPEN_ACCESS:Y AND FIRST_PDATE:[2020-01-01 TO 2026-12-31]` (9,062 hits) and the first 1,200 records were sampled; 992 full texts could be retrieved. Sentences containing "palmitoyl" and a pattern of the form `Cys<number>` or `C<number>` were extracted, with the site pattern required not to be preceded by a letter or digit (so that `zDHHC1` is not read as site 1). A protein token was taken from the 40 characters to the left or right of the match and kept only if it looked like a gene symbol (contains a digit, is upper case, or matches a lower-case prefix followed by an upper-case letter). Tokens were resolved against a local dictionary of gene names from the 28,711 reviewed sequences already downloaded; a token was accepted only when it mapped to a single accession (exact match, or a unique prefix). A mention became a positive only if the residue at that position in the canonical sequence is a cysteine. Negatives were up to six other cysteines of the same protein that are neither mentioned in the set nor annotated as S-palmitoylated in UniProt 2026_03. A position-specific log-odds motif trained on the UniProt data was used as an independent check of whether the extracted positives carry a palmitoylation motif.

### 4.11 Proteome-wide scan and decoy construction

All 262,286 cysteines of the reviewed human sequences were scored in chunks with the forest view of the released model, calibrated on the same out-of-fold validation rows. The ESM view was not run at this scale. The null was 20,000 shuffled-flank decoys that keep the central cysteine but destroy the motif; a C→A decoy would remain motif-shaped. FDR at a threshold was estimated as (decoys at or above the threshold) × (targets / decoys) divided by (targets at or above the threshold), with q-values as the running minimum over lower thresholds. Because the decoy tail is thin and the shuffled null is not proven to represent unmodified cysteines, the tiers are candidate rankings rather than demonstrated error rates.

### 4.12 AlphaFold features and paired evaluation

Human proteins with at least one annotated site and a random sample of background human proteins were queried through the AlphaFold DB prediction API; models were retrieved for chains shorter than 700 residues, covering 1,950 human sites (463 positives in 476 proteins). Per site we recorded the C-alpha pLDDT, mean pLDDT within ±5 positions, and C-alpha atoms within 10 Å excluding neighbours closer than four positions. The same random forest and cluster-disjoint folds were used with and without the features, so the comparison is paired at the fold level. The prediction service exposes the sequence path only.

### 4.13 Reproducibility

The pipeline is organised as 25 scripts in `src/` (pipeline table in the repository README, ordered commands in `docs/REPRODUCE.md`) with fixed seeds; processed data, split assignments, out-of-fold predictions and per-fold metrics are included. Runtime on a 16-core CPU machine is approximately 3 h end-to-end, of which ESM embedding extraction is the largest component.

---

## Data and code availability

The corrected v2 code, processed data (including the canonical fold assignments), out-of-fold predictions and proteome-wide scores are available at https://github.com/opus940810/palm-dl (release tag v2.1; the exact commit, file hashes and reproduction commands are listed in docs/REPRODUCE.md). The earlier v1.0.0 release is archived at Zenodo (DOI 10.5281/zenodo.23455873); it predates the label correction and must not be used for the results reported here, and a corrected Zenodo version will be deposited before publication. The prediction service is available at [SERVICE-URL]. Supplementary Tables S1–S7 accompany the manuscript. Raw data are from UniProt release 2026_03 and are freely available under CC BY 4.0; AlphaFold models are from the AlphaFold DB (CC BY 4.0).

## Competing interests

The authors declare no competing interests.

## Author contributions

[To be completed.]

## Acknowledgements

[To be completed.]

---

## References

1. Mesquita FS, Abrami L, Linder ME, Bamji SX, Dickinson BC, van der Goot FG. Mechanisms and functions of protein S-acylation. *Nat Rev Mol Cell Biol.* 2024;25:488–509. doi:10.1038/s41580-024-00700-8
2. Anwar MU, van der Goot FG. Refining S-acylation: structure, regulation, dynamics, and therapeutic implications. *J Cell Biol.* 2023;222:e202307103. doi:10.1083/jcb.202307103
3. Shang J, Ding M, Zhou X. Recent advances in S-palmitoylation and its emerging roles in human diseases. *J Hematol Oncol.* 2025;18:83. doi:10.1186/s13045-025-01738-7
4. Zhou F, Xue Y, Yao X, Xu Y. CSS-Palm: palmitoylation site prediction with a clustering and scoring strategy (CSS). *Bioinformatics.* 2006;22(7):894–896. doi:10.1093/bioinformatics/btl013
5. Ren J, Wen L, Gao X, Jin C, Xue Y, Yao X. CSS-Palm 2.0: an updated software for palmitoylation sites prediction. *Protein Eng Des Sel.* 2008;21(11):639–644. doi:10.1093/protein/gzn039
6. Xue Y, Chen H, Jin C, Sun Z, Yao X. NBA-Palm: prediction of palmitoylation site implemented in Naïve Bayes algorithm. *BMC Bioinformatics.* 2006;7:458. doi:10.1186/1471-2105-7-458
7. Weng SL, Kao HJ, Huang CH, Lee TY. MDD-Palm: identification of protein S-palmitoylation sites with substrate motifs based on maximal dependence decomposition. *PLoS One.* 2017;12(6):e0179529. doi:10.1371/journal.pone.0179529
8. Ning W, Jiang P, Guo Y, Wang C, Tan X, Zhang W, et al. GPS-Palm: a deep learning-based graphic presentation system for the prediction of S-palmitoylation sites in proteins. *Brief Bioinform.* 2021;22(2):1836–1847. doi:10.1093/bib/bbaa038
9. Li S, Yu K, Wu G, Zhang Q, Wang P, Zheng J. pCysMod: prediction of multiple cysteine modifications based on deep learning framework. *Front Cell Dev Biol.* 2021;9:617366. doi:10.3389/fcell.2021.617366
10. Bandyopadhyay SS, et al. RFCM-PALM: in-silico prediction of S-palmitoylation sites in the synaptic proteins for male/female mouse data. *Int J Mol Sci.* 2021;22(18):9901. doi:10.3390/ijms22189901
11. Blanc M, David F, Abrami L, Migliozzi D, Armand F, Bürgi J, et al. SwissPalm: protein palmitoylation database. *F1000Res.* 2015;4:261. doi:10.12688/f1000research.6464.1
12. Meng Y, Zhang L, Zhang L, Wang Z, Wang X, Li C. CysModDB: a comprehensive platform with the integration of manually curated resources and analysis tools to facilitate the study of cysteine modifications. *Brief Bioinform.* 2022;23(6):bbac460. doi:10.1093/bib/bbac460
13. Meng L, Chan W, Huang L, Liu L, Chen X, Zhang W. Mini-review: recent advances in post-translational modification site prediction based on deep learning. *Comput Struct Biotechnol J.* 2022;20:3522–3532. doi:10.1016/j.csbj.2022.06.045
14. Deng M, et al. Deep-Palm: an integrated deep learning framework for structure-aware prediction of protein S-palmitoylation. *bioRxiv* 2026. doi:10.64898/2026.03.05.709753 (preprint)
15. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res.* 2025;53:D609–D617. doi:10.1093/nar/gkae1010

---

## Figure legends

**Figure 1. Dataset composition.** (a) Organisms contributing the largest numbers of annotated S-palmitoyl cysteine sites (human and mouse highlighted; remaining 784 organisms pooled). (b) Evidence tier of positive sites (log scale): 455 experimentally supported versus 3,969 inferred, mostly by similarity. (c) Annotated topology of positive sites.

**Figure 2. Split protocol determines measured performance.** A single random forest was trained and evaluated under four protocols. (a) AUC and class-balanced MCC. (b) Percentage of test positives whose exact 31-residue window is present among training positives. Random site-level and protein-level splits leak 41–45% of windows; cluster and species splits leak none, for positives and negatives alike.

**Figure 3. Discrimination under cluster-held-out cross-validation.** (a) ROC and (b) precision–recall curves pooled over out-of-fold predictions; the table reports unweighted fold means, so the pooled curves and the table are related but not identical estimands. Models: PWM, logistic regression, random forest, CNN, ESM-2 view and nested PalmDL (4,424 positives, 67,266 sites).

**Figure 4. Performance on human sites under three evaluation regimes.** (a) AUC and (b) balanced MCC for multi-species training with human test folds, human-only training, and the strict species holdout. The species-holdout bars are single-model (random forest/CNN) results and are not a full ensemble; negative compositions differ between the three regimes.

**Figure 5. Calibration.** Reliability diagrams for the random forest, CNN, ESM-2 view and calibrated ensemble at natural prevalence.

**Figure 6. Sequence motif learned from the data.** Position-specific log-odds of amino acids from −7 to +8 around the central cysteine, positives versus background negatives (log2 scale).

**Figure 7. Cross-species transfer is heterogeneous.** AUC of single-species PalmDL models, rows = training species, columns = test species, for the six species with at least 85 S-palmitoyl cysteine sites. Most mammalian pairs transfer at 0.83–0.95 with a pig→human exception; the plant model reaches 0.74–0.78 on mammals. Empty cells have fewer than 20 eligible test positives after cluster removal. Diagonal entries come from a protein-cluster-disjoint 20% holdout of the same species.

**Figure 8. What the model uses and how far its predictions can be trusted.** (a) Mean change in predicted probability when each position of 300 human positive windows is mutated to each amino acid (forest view). Cysteine at ±1 increases the score (cysteine clusters), acidic residues at −1/−2 decrease it, and the central position is neutral because every candidate is a cysteine. (b) Validation-target operating points: precision against recall on the outer test folds when the threshold is chosen on validation at empirical targets of 1%, 5%, 10% and 20%; observed held-out false-discovery proportions are 4.7%, 10.0%, 15.1% and 24.1%.

**Figure 9. Proteome-wide scan and candidate shortlist.** (a) Score distributions of all 262,286 human cysteines and of 20,000 shuffled-flank decoys (log scale). (b) Estimated FDR as a function of shortlist size, with the five reported tiers marked. (c) Composition of each tier: known UniProt sites versus novel cysteines.

## Table legends

**Table 1.** Cluster-held-out cross-validation (all species, 5 folds, 67,266 sites). AUC, AUPR, MCC and balanced MCC with 95% bootstrap CIs for PalmDL.

**Table 2.** Leakage audit: identical-window overlap and performance under four splitting protocols.

**Table 3.** Human evaluations: multi-species training, human-only training, species holdout, and the hard-negative benchmark.

**Table 4.** Comparison with published S-palmitoylation predictors.

**Table 5.** Cross-species transfer matrix (AUC) for six species.

**Table 6.** Error structure of the held-out ensemble predictions by evidence tier and annotated topology.

**Table 7.** Automated literature mining as an external benchmark: yields, motif consistency and model performance under progressively stricter extraction filters.

**Table 8.** Proteome-wide scan of the reviewed human proteome at five candidate tiers, with recovery of known sites and the fraction of novel predictions.

**Table 9.** AlphaFold structural features on human sites: paired cross-validated comparison and feature distributions.
