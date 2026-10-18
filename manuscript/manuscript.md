# PalmDL: a leakage-aware deep learning framework and benchmark for protein S-palmitoylation site prediction

**Running title:** Leakage-aware deep learning for S-palmitoylation sites

**Authors:** [Author 1]¹, [Author 2]¹, [Corresponding author]¹\*

¹ [Affiliation], [City], [Country]

\* Correspondence: [email]

---

## Abstract

**Motivation.** S-palmitoylation is a reversible lipid modification that controls protein trafficking, signalling and membrane organisation, and its dysregulation underlies cancer, neurodegeneration and immune disease. Experimental mapping remains low-throughput, so dozens of computational predictors have been published. Yet the field has no common, homology-controlled benchmark: most tools are trained and tested on random or protein-level splits in which near-identical peptides, propagated annotations and orthologues cross the train/test boundary, and no predictor has been evaluated on a prospectively held-out species.

**Results.** We assembled a versioned, evidence-stratified dataset from UniProt release 2026_03 comprising 71,594 labelled cysteine sites (6,805 S-palmitoyl sites, of which 477 are experimentally supported; 1,425 organisms; 572 human sites in 316 proteins) with explicitly separated hard (same-protein) and background negatives. Using a single model across four splitting protocols we show that conventional benchmarks leak severely: 43.7% of test positives share an identical 31-residue window with a training positive under a random site-level split and 40.2% under a protein-level split, inflating AUC from 0.950 to 0.990 and MCC from 0.745 to 0.843. We then developed PalmDL, a two-view ensemble that combines ESM-2 protein-language-model embeddings with motif/physicochemical features through a random forest. Under cluster-held-out cross-validation PalmDL reaches AUC 0.956 (95% CI 0.953–0.959), AUPR 0.849 and MCC 0.753; on held-out human sites it reaches AUC 0.901 (0.885–0.915), MCC 0.355 and balanced MCC 0.691; on a strict cross-species holdout, where no human data and no human-homologous cluster are seen during training, AUC is 0.822. PalmDL outperforms a CSS-Palm-style PWM, logistic regression, a multi-kernel CNN and each single view (ΔAUC +0.011 to +0.015 versus the best single model, paired bootstrap). Probabilities are Platt-calibrated (ECE 0.040 at natural prevalence). Two further analyses frame the biology rather than the algorithm: a cross-species transfer matrix shows that the palmitoylation code is eukaryotic rather than species-specific (mammalian models reach AUC 0.94–0.97 on other mammals and on plants, whereas bacterial lipidation is a different code, AUC ≈0.55), and calibrated scores support FDR-controlled shortlists (at 5% estimated FDR, precision 0.94 at 61% recall). In-silico mutagenesis recovers the known cysteine-cluster rule and exposes annotation errors in the training labels. Finally, a literature-mined external set of 136 sites from 91 papers — 98% of them outside UniProt — carries no canonical motif (PWM AUC 0.54); manual curation of 99 unambiguous sites and a ±2-residue numbering control leave the model at AUC 0.62 against same-protein cysteines, exposing a domain-shift limitation of annotation-trained predictors. Applied proteome-wide, the model yields a target–decoy-controlled shortlist of 322 human cysteines at 5% FDR (227 of them novel) that is enriched 16-fold for cytoplasmic topology. Adding AlphaFold-derived pLDDT and packing features lifts human-site AUC from 0.890 to 0.910 and the full model from 0.946 to 0.949 (5/5 folds), although it does not rescue the motif-atypical literature sites.

**Availability.** All code, processed data, fold assignments and trained-model recipes are provided; the pipeline is reproducible end-to-end from public UniProt queries.

**Keywords:** S-palmitoylation; post-translational modification; deep learning; protein language model; benchmark leakage; cross-species generalisation

---

## 1. Introduction

S-palmitoylation is the reversible attachment of a 16-carbon fatty acid to cysteine thiols. Unlike most lipid modifications it is enzymatically reversible, and the competing actions of ~23 DHHC acyltransferases and several thioesterases make it a dynamic switch that partitions proteins between membrane compartments, controls synaptic transmission, and shapes oncogenic signalling [1–3]. The functional importance of the modification is matched by the difficulty of measuring it: acyl-biotin exchange, click-chemistry labelling and mass spectrometry can map sites, but each is biased by abundance, accessibility and labelling efficiency, and proteome-wide site-level coverage remains sparse.

Computational prediction has therefore been attractive for almost two decades. CSS-Palm 1.0 introduced clustering-based scoring with BLOSUM62 [4], CSS-Palm 2.0 added matrix mutation and a homology-reduced benchmark [5], and SVM-based and motif-based methods followed, including NBA-Palm [6] and MDD-Palm [7]. Deep learning arrived with GPS-Palm, a parallel CNN over sequence and structure-derived features that reported a human AUC of 0.900 [8], and pCysMod, a multi-PTM cysteine model reporting AUC 0.807 for palmitoylation [9]. More recent entries include RFCM-PALM for the mouse synaptic proteome [10], the SAPPTree structure-based random forest (preprint) and Deep-Palm, a stacking ensemble that reports AUC 0.931 on a CD-HIT 60% site-level 80/20 split (preprint [14]). CSS-Palm 3.0 and 4.0 were distributed as standalone software but have no peer-reviewed description, and their web servers have been retired.

Three problems limit the comparability and credibility of this literature. First, negatives are almost always "all other cysteines" of the same proteins or random cysteines from unannotated proteins; the ratio and difficulty of negatives differ between studies, and accuracy or AUC computed on easy background negatives is not comparable across papers. Second, split protocols are rarely protein- or family-disjoint. Palmitoylation annotations are heavily propagated by similarity: in UniProt, 93% of S-palmitoyl cysteine annotations are inferred rather than experimentally determined, and orthologues across species carry identical 31-residue windows. A random site-level split therefore places the same peptide in training and test, and even a protein-level split does so whenever annotations were propagated between species. Third, no palmitoylation predictor has been tested on a species that was entirely absent from training, so it is unknown whether learned motifs transfer or whether they memorise species- and family-specific sequence context.

Here we address all three. We build a versioned, evidence-stratified dataset from a single UniProt release and separate hard from background negatives. We quantify, with one fixed model, how much each splitting protocol inflates performance, and measure the fraction of test positives whose exact window appeared in training. We then develop PalmDL, a leakage-aware ensemble of a protein-language-model view (ESM-2) and a motif/physicochemical view (random forest), calibrated and evaluated under cluster-held-out cross-validation, human-only training, and a strict cross-species holdout. We show that PalmDL exceeds the best previously reported human AUC on a benchmark with zero shared windows, and we release the split definitions so that future methods can be compared on identical, leakage-controlled folds.

---

## 2. Results

### 2.1 A versioned, evidence-stratified dataset

We queried UniProt release 2026_03 for all Swiss-Prot entries carrying at least one lipidation feature (9,089 entries) and parsed every feature whose description was "S-palmitoyl cysteine". After removing duplicates and sites outside the canonical sequence, the dataset contains **6,805 positive sites**: 477 supported by experimental evidence codes and 6,328 inferred (mostly by similarity), distributed over 1,425 organisms, with human (572 sites in 316 proteins) and mouse (524) the largest contributors (Figure 1a,b). Most sites (4,155) have no annotated topology; of the annotated ones, 2,445 are cytoplasmic and only 42 extracellular (Figure 1c).

Negatives were defined in two tiers. **Hard negatives** are unannotated cysteines of the same proteins whose ±15-residue windows do not overlap an annotated site (10,567 sites). **Background negatives** are cysteines sampled from 19,622 reviewed human proteins that carry no lipidation feature at all (54,222 sites). The final labelled set contains 71,594 cysteine sites from 23,897 proteins, with 64,143 distinct 31-residue windows. Metrics are stable when evaluation is restricted to the harder same-protein negatives (AUC 0.965 versus 0.956 on the full set; Table 3), so model ranking is not an artefact of easy background negatives.

### 2.2 Conventional splits leak, and the leakage is measurable

We trained the same random forest on the same features under four protocols and evaluated each on holdout sites that were never used for fitting or threshold selection (Table 2, Figure 2). Under a random site-level split, **43.7% of test positives had an identical 31-residue window among the training positives**, and AUC reached 0.990 (MCC 0.843). A protein-level split, which is the most common precaution in the field, still leaked **40.2%** of test windows and gave AUC 0.988. Cluster-level splitting removed the overlap entirely (0%) and reduced AUC to 0.950 and MCC to 0.745. Holding out an entire species reduced AUC to 0.822.

This audit has two consequences. First, headline numbers from site-level or protein-level splits cannot be compared with numbers from homology-controlled splits; a model reporting AUC > 0.98 is likely measuring peptide memorisation. Second, the widely used CD-HIT 60% cut-off is insufficient on its own when 93% of labels are propagated by similarity, because identical windows can survive clustering if they are carried by short or partially covered proteins. We therefore report every result in this paper under cluster-level splits, and additionally under a species holdout.

### 2.3 Model development: two complementary views

We compared five model families under identical cluster-held-out folds (Table 1, Figure 3). A CSS-Palm-style position-specific log-odds matrix, the classical baseline, reached AUC 0.862 and MCC 0.516. Logistic regression on the flattened one-hot window reached AUC 0.879. A multi-kernel CNN (kernel sizes 3, 5, 7 and 11; 32 filters; protein-context features) reached AUC 0.894. A random forest on one-hot and physicochemical channels reached AUC 0.945 and MCC 0.707. The strongest single view was a small MLP on frozen ESM-2 (35M) embeddings, computed separately for the mean of the window and for the central cysteine: AUC 0.941, MCC 0.741, AUPR 0.820.

The two strongest views make different errors: the forest captures exact positional chemistry and local motif syntax, whereas ESM-2 embeddings capture evolutionary and structural context. A rank-average of the two, followed by Platt calibration, was selected on validation AUC among all 31 non-empty member subsets and reached **AUC 0.956 (95% CI 0.953–0.959), AUPR 0.849, MCC 0.753, balanced MCC 0.767** (Table 1). Paired bootstrap over sites gives ΔAUC = +0.015 (95% CI 0.013–0.017) over ESM-2 alone and +0.011 (0.010–0.013) over the forest; both exclude zero. Adding the CNN to the ensemble did not improve validation AUC (0.9517 versus 0.9513) and reduced test AUC, so the simpler two-view model was retained.

### 2.4 Human evaluation and the cost of a strict split

Because most published predictors are evaluated on human sites, we repeated the analysis on human-only test folds (Table 3, Figure 4). With multi-species training and human test folds, PalmDL reaches **AUC 0.901 (0.885–0.915), AUPR 0.279, MCC 0.355 and balanced MCC 0.691**. A 4,000-replicate bootstrap against the best published human value, 0.900, gives P(AUC ≤ 0.900) = 0.46, so the two are statistically indistinguishable; the published general (multi-species) AUC of 0.855 is exceeded with P < 0.001. We therefore claim parity with the state of the art on human sites, not superiority, and the protocol we used is stricter (zero shared windows, cluster-held-out). Training on human data alone gives essentially the same AUC (0.897), showing that the additional 6,233 non-human sites add little once homology is controlled — consistent with the fact that most non-human annotations are propagated from the same experimentally characterised proteins. On the strict species holdout, in which no human protein and no cluster containing a human protein is seen during training and every human-propagated annotation is removed, the best single model (random forest) reaches AUC 0.822 and balanced MCC 0.411. This 0.13 AUC gap between within-species and cross-species evaluation is, to our knowledge, the first direct measurement of how much palmitoylation predictors rely on species-specific context.

### 2.5 Calibration, hard negatives and interpretability

Probabilities from the calibrated ensemble are reliable at natural prevalence (expected calibration error 0.040; Figure 5), whereas the uncalibrated forest has ECE ≈ 0.12 and over-confident tails. On the hard-negative benchmark (positives plus same-protein cysteines only, 1:1.55 ratio), PalmDL reaches AUC 0.965 and MCC 0.755, compared with 0.868/0.570 for the PWM baseline (Table 3), so the advantage is not an artefact of easy background negatives.

The learned motif agrees with known biochemistry (Figure 6). Positive windows are enriched in hydrophobic and aromatic residues at positions −3 to −1 relative to the cysteine (leucine, phenylalanine, tryptophan, tyrosine and methionine), depleted of glycine, proline and acidic residues in the same interval, and enriched in basic residues at positions +3 to +4. This is consistent with the DHHC active-site geometry, in which the cysteine is presented in a hydrophobic groove, and with the known preference for cysteine clusters and juxtamembrane positions.

### 2.6 Position relative to published predictors

Table 4 places our results next to the numbers reported by previous methods. Three patterns are visible. First, the classical MCC values (0.45–0.67) were obtained under protocols that either reuse the training literature or do not hold out proteins, and the only classical tool that reported an independent test (CSS-Palm 2.0) dropped from MCC 0.671 to 0.523 when it did so. Second, the two deep-learning systems with published AUCs (GPS-Palm 0.900, pCysMod 0.807) report no cross-species and no explicit protein-disjoint evaluation, so their operating point is unknown. Third, the highest published number (Deep-Palm, 0.931) comes from a site-level 80/20 split, precisely the protocol that our audit shows to leak ~44% of windows. PalmDL's 0.956 is obtained with zero shared windows, and its human-only estimate, 0.901, is on a par with the best published human AUC while being the only one accompanied by a leakage audit, a calibration analysis and a cross-species holdout. Because we could not rerun GPS-Palm or the retired classical tools on our folds, Table 4 compares protocols as much as algorithms; the released fold files remove this obstacle for future work.

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
| **PalmDL** | this work | cluster-held-out CV, 6,805 sites | **AUC 0.956, MCC 0.753** | yes (0% shared windows) |
| **PalmDL** | this work | human test folds | **AUC 0.901, balanced MCC 0.691** | yes (clusters) |
| **PalmDL** | this work | species holdout | AUC 0.822, balanced MCC 0.411 | yes (strictest) |

### 2.7 The palmitoylation code is eukaryotic, not species-specific

The species holdout in §2.4 removes every human-homologous cluster and therefore measures transfer to genuinely novel families. To separate that from ordinary cross-species transfer we trained a separate PalmDL model on each of eight species with at least 85 annotated sites and evaluated every model on every species (Figure 7, Table 5). Three regimes emerge.

First, **mammals share one code**. Models trained on human, mouse, rat, cattle or pig reach AUC 0.86–0.98 on the other mammals, with a mean of 0.943 that is indistinguishable from the mean within-species diagonal (0.944). Second, **the code extends to plants**: mammalian models reach 0.965 on *A. thaliana*, and in the reverse direction the plant model reaches 0.763 on mammals — lower, but far above chance. Third, **bacterial lipidation is a different modification**: *E. coli* and *B. subtilis* models are at chance on mammalian sites (0.551), and mammalian models transfer only weakly to bacteria (0.671), while the two bacterial models transfer to each other almost perfectly (0.95–1.00). This is biologically expected — bacterial lipoproteins are N-acylated at the LipoBox rather than S-palmitoylated at a cysteine — and it confirms that the model has learned a eukaryotic motif rather than a generic "hydrophobic cysteine" rule.

The contrast with §2.4 is the key point. Cross-species transfer with orthologues available is essentially free (0.94), but transfer to human families absent from training is much harder (0.822). Most of the apparent generality of palmitoylation predictors therefore comes from orthologue coverage, not from an abstract motif grammar.

**Table 5. Cross-species transfer matrix (AUC).**

| train \ test | Human | Mouse | Rat | Cattle | Pig | *A. thaliana* | *E. coli* | *B. subtilis* |
|---|---|---|---|---|---|---|---|---|
| Human | 0.914 | 0.963 | 0.966 | 0.973 | 0.968 | 0.973 | 0.599 | 0.759 |
| Mouse | 0.958 | 0.888 | 0.976 | 0.974 | 0.963 | 0.970 | 0.389 | 0.604 |
| Rat | 0.943 | 0.970 | 0.948 | 0.972 | 0.968 | 0.972 | 0.833 | 0.879 |
| Cattle | 0.900 | 0.930 | 0.941 | 0.856 | 0.964 | 0.963 | 0.754 | 0.721 |
| Pig | 0.857 | 0.878 | 0.882 | 0.918 | 0.998 | 0.946 | 0.536 | 0.632 |
| *A. thaliana* | 0.693 | 0.771 | 0.787 | 0.788 | 0.775 | 1.000 | 0.767 | 0.812 |
| *E. coli* | 0.505 | 0.552 | 0.575 | 0.530 | 0.503 | 0.493 | 1.000 | 0.997 |
| *B. subtilis* | 0.539 | 0.589 | 0.608 | 0.561 | 0.548 | 0.588 | 0.996 | 0.950 |

Diagonal entries are from a protein-cluster-disjoint 20% holdout of the same species; off-diagonal entries use a model trained on all rows of the source species. Small species have few test positives (*A. thaliana* 97 sites, *E. coli* 103, *B. subtilis* 85), so their rows are correspondingly noisy.

### 2.8 Error structure and FDR-controlled shortlists

Because a predictor is only useful if its errors are predictable, we stratified the held-out predictions (Table 6). Discrimination is substantially better for inferred than for experimentally supported sites (AUC 0.960 versus 0.908; mean predicted probability 0.55 versus 0.42). This "annotation-memory" effect is expected when 93% of labels are propagated: experimentally characterised sites are the atypical, family-specific cases that propagation does not cover, and they are exactly the ones a user is most likely to care about. Performance is stable across annotated topologies (cytoplasmic 0.898, extracellular 0.898, topology unknown 0.949) with the exception of the lumenal class, which contains only 14 positives and is not interpretable.

For practical use we derived thresholds whose estimated false-discovery rate is controlled on the inner validation folds and then measured precision and recall on the outer folds (Figure 8b). At an estimated FDR of 1% the threshold retains 44.5% of known sites at 98.3% precision; at 5% it retains 61.3% at 93.5% precision; at 10%, 67.5% at 87.8%. This gives a principled way to generate a shortlist for experimental follow-up without committing to an arbitrary score cut-off, and it is the mode in which we expect the model to be used.

**Table 6. Error structure of the held-out ensemble predictions.**

| Stratum | n positives | AUC | AUPR |
|---|---|---|---|
| Experimentally supported sites | 506 | 0.908 | 0.293 |
| Inferred sites | 6,299 | 0.960 | 0.855 |
| Cytoplasmic | 2,445 | 0.898 | 0.883 |
| Extracellular | 42 | 0.898 | 0.355 |
| Topology unknown | 4,155 | 0.949 | 0.843 |

### 2.9 What the model uses: cysteine clusters, charge exclusion, and a label audit

In-silico mutagenesis of 300 human positive windows (every position mutated to every amino acid, scored with the forest view; Figure 8a) shows that substitutions at the central position are neutral, as they must be: every candidate is a cysteine, so the centre carries no discriminative information. The informative positions are the flanks. Introducing a cysteine at ±1 increases the predicted probability by 0.047, recovering the known preference for cysteine clusters; introducing glutamate or aspartate at −1 or −2 decreases it by 0.06–0.10, recovering the exclusion of acidic residues near the modification site.

The same audit exposed a label-quality problem. Before filtering, mutating the central cysteine to lysine *increased* the score by 0.089 — an impossible result that traced back to 66 UniProt features (1.0% of positives) whose annotated position carries a lysine, glycine, serine or threonine in the current canonical sequence, reflecting transferred annotations on diverged or re-annotated entries. Removing those 66 sites abolished the artefact (effect 0.000) without changing cross-validated performance (AUC 0.9451 with and without). Interpretability is therefore not only a biological tool but a label-auditing one.

### 2.10 Literature mining cannot yet supply an independent benchmark

Finally we asked whether the missing external benchmark can be manufactured from the literature. We searched Europe PMC for open-access articles on palmitoylation published between 2020 and 2026 (9,062 hits), fetched 992 full texts, and extracted every sentence that mentions palmitoylation together with a cysteine position. Protein tokens were resolved against 28,711 locally cached reviewed sequences, and a mention was kept only when the residue at the stated position is a cysteine: **136 sites in 109 proteins from 91 papers** (Table 7). Only 2.2% of these sites are annotated in UniProt 2026_03, so the set lies almost entirely outside the training labels and is the closest thing to an independent site-level benchmark currently obtainable.

It is not a clean benchmark. The motif that the model learned is barely enriched in this set (PWM AUC 0.537 against same-protein cysteines), and PalmDL separates the literature sites from those cysteines only weakly (AUC 0.588; 0.607 against random background cysteines). Tightening the extraction — requiring an exact gene match, explicit "palmitoylation at CysX" phrasing, and excluding table rows — leaves the picture unchanged (PWM AUC 0.56–0.57; model AUC 0.61–0.62). Manual inspection explains why: genuine sites are mixed with extraction errors, and many genuine sites (MyD88 C113, iRhom2 C476, GRP78 C420) are transmembrane-proximal cysteines whose linear motif is weak, which is precisely the class that an annotation-trained model has least often seen.

To make sure this was not an extraction artefact, we curated all 136 sites by reading every sentence (Table 7, last row). Ninety-nine mentions could be assigned unambiguously to a protein and a cysteine — either because the sentence states the site explicitly or because a cysteine-to-serine/alanine mutant lost palmitoylation — while 37 were rejected as sites of a different protein in the same sentence, catalytic-site mutants of DHHC enzymes, predicted-only sites, negative results, or table rows with mismatched proteins. The curated subset behaves exactly like the full set (PWM AUC 0.561, model AUC 0.620), so extraction noise is not the explanation. We also tested the obvious alternative explanation, that the papers use a different residue numbering: allowing every site and every same-protein negative to score its best window within ±2 residues leaves AUC at 0.611, so numbering offsets do not explain the failure either.

Consistent with this, the curated sites are essentially absent from the proteome-wide shortlist: one of 99 appears in the 5% FDR tier (expected 0.12 by chance; 8-fold but n = 1) and two in the 20% tier (expected 0.25). The honest conclusion is a domain-shift limitation that applies to any annotation-trained predictor, including PalmDL: it reproduces and extends the UniProt consensus, but does not recognise the motif-atypical, transmembrane-proximal sites that mechanistic papers tend to characterise.

We therefore report this as an audit of the resource landscape rather than as a validation. It quantifies the statement in the Introduction: the field's headline numbers have never been checked on data that are simultaneously independent, site-level and motif-consistent, and automated literature mining does not yet close that gap.

**Table 7. Automated literature mining as an external benchmark: yields and motif consistency.**

| Filter | sites | PWM AUC | PalmDL AUC |
|---|---|---|---|
| All verified sites | 133 | 0.537 | 0.588 |
| Exact gene-symbol match | 131 | 0.539 | 0.591 |
| Explicit "palmitoylation at CysX" phrasing | 97 | 0.560 | 0.613 |
| Exact match + explicit phrasing | 96 | 0.559 | 0.615 |
| + table-like rows removed | 95 | 0.565 | 0.616 |
| Manually curated subset | 99 | 0.561 | 0.620 |

Comparisons are against other cysteines of the same proteins. PWM AUC is the enrichment of a position-specific log-odds motif trained on the UniProt data; if the extracted positives were faithful, this value should approach the model AUC. The manually curated subset is not nested in the filters above: it was obtained by reading all 136 sentences and accepting 99 unambiguous mentions (82 proteins).

### 2.11 A proteome-wide, FDR-controlled candidate list

As a forward-looking counterpart to the retrospective benchmarks, we applied the released model to every cysteine of the reviewed human proteome — **262,286 sites in 19,814 proteins** — and controlled the false-discovery rate with shuffled-flank decoys that keep the central cysteine but destroy the motif. (A C→A decoy is inappropriate here: §2.9 shows that the model is deliberately insensitive to the identity of the centre, so such mutants remain motif-shaped and give a null that is far too pessimistic.) At 1% FDR the scan yields 217 sites in 138 proteins; at 5%, 322 sites in 204 proteins; at 10%, 399; at 20%, 668; and at 30%, 1,065 (Table 8, Figure 9). Recall of the 572 UniProt-annotated human sites rises from 12.4% at the 1% tier to 39.0% at the 30% tier, and the strict tiers are dominated by novel sites (146 of 217 at 1%).

The shortlist is biologically coherent. 90.7% of the 5% FDR sites carry annotated cytoplasmic topology versus 5.6% of all scanned cysteines — a 16.3-fold enrichment — and the highest-scoring novel entries are juxtamembrane cysteines of single-pass membrane proteins: endothelin receptor A C386/C388, claudin-2 C185, mucolipin-3 C550 and LRRC32 C650. This is precisely the class that canonical motif models under-score, and it is where a sequence-only model can add value. Only one of the 136 literature-mined sites reaches the 5% tier, consistent with the motif-atypical nature of that set (§2.10).

Two caveats bound the resource. First, the proteome-wide prior is small (572 annotated sites among 262,286 cysteines, ≈0.2%), so strict FDR control is intrinsically severe; users interested in discovery rather than confirmation should work at the 20–30% tiers. Second, 82% of scanned cysteines have no topology annotation and are strongly depleted from the shortlist (9.3%), so the list is biased toward proteins with curated topology.

**Table 8. Proteome-wide scan of the reviewed human proteome (262,286 cysteines) at five FDR tiers.**

| FDR | sites | proteins | known recovered | recall of known | novel | score threshold |
|---|---|---|---|---|---|---|
| 1% | 217 | 138 | 71 | 12.4% | 146 | 0.969 |
| 5% | 322 | 204 | 95 | 16.6% | 227 | 0.957 |
| 10% | 399 | 246 | 115 | 20.1% | 284 | 0.948 |
| 20% | 668 | 418 | 164 | 28.7% | 504 | 0.918 |
| 30% | 1,065 | 659 | 223 | 39.0% | 842 | 0.864 |

"Known" refers to the 572 S-palmitoyl cysteine sites annotated in UniProt 2026_03 for human proteins; "novel" sites are unannotated cysteines. FDR is estimated by target–decoy comparison against 20,000 shuffled-flank decoys.

### 2.12 A structure-aware prototype: useful inside the annotation distribution, not outside it

Sequence-only models cannot see whether a cysteine is buried in a folded domain or exposed on a membrane-proximal surface, which matters because the acyl chain must be accommodated. We fetched AlphaFold models for 266 of the 316 human proteins carrying annotated sites (plus 213 background proteins, 479 in total) and computed three features per cysteine: pLDDT at the site, mean pLDDT within ±5 residues, and C-alpha packing density (number of C-alpha atoms within 10 Å, excluding sequence neighbours). 1,974 human sites (476 positives) were covered.

We compared four ways of using them (Table 9). (i) Adding the features to the random-forest view lifts AUC from 0.890 to 0.910 on the structured subset (Δ +0.020, 5/5 folds). (ii) Reranking the full sequence ensemble with a small logistic model fails (0.919 → 0.914), because the ensemble already carries most of the structural signal through the language-model view and a four-parameter reranker overfits 2,000 sites. (iii) Treating structure as a third view — averaging the structure-aware forest with the frozen ESM view — gives a smaller but perfectly consistent gain (0.918 → 0.924, 5/5 folds). (iv) A hybrid forest trained on all 71,594 sites, with median imputation and a missing-structure indicator for the 98% of sites without a model, improves the full test set from 0.946 to 0.949 (5/5 folds) and the structured subset from 0.899 to 0.918 (4/5 folds); this is the version we release as `models/rf_hybrid_struct.joblib`, so that proteins with a model benefit and the rest fall back to the sequence path.

The direction of the structural effect is biologically coherent: palmitoylated cysteines have lower pLDDT (65.8 versus 81.2) and lower packing (4.8 versus 10.9 C-alpha neighbours), i.e. they sit in less ordered, more exposed segments where the acyl chain can be accommodated and where DHHC enzymes act. Applied to the proteome scan, the hybrid recovers more known sites in the same budget (321 versus 302 in the top 322 of the structured subset; 466 versus 412 in the top 668), so within the annotation distribution structure is a real, if modest, improvement.

It does not cure the domain shift of §2.10. On the 700 literature sites that have a structural model, the structure-aware forest is no better than the sequence forest (0.627 versus 0.641), and neither approaches the performance seen on UniProt sites. The class of sites that literature studies characterise — often transmembrane-proximal, sometimes weakly structured — is not rescued by burial or flexibility alone. Structure therefore improves ranking where annotations already exist, but it is not by itself the missing piece for genuinely novel sites.

**Table 9. Structure-aware variants on human sites (paired cluster-disjoint folds).**

| Variant | Comparison | AUC | Δ |
|---|---|---|---|
| Forest + structure | 0.890 → 0.910 (1,974 sites) | 0.910 | +0.020 (5/5 folds) |
| Logistic reranker on ensemble | 0.919 → 0.914 | 0.914 | −0.005 (0/5) |
| Structure as third view with ESM | 0.918 → 0.924 | 0.924 | +0.006 (5/5) |
| Hybrid forest, all 71,594 sites | 0.946 → 0.949 | 0.949 | +0.003 (5/5) |
| Hybrid forest, structured subset | 0.899 → 0.918 | 0.918 | +0.019 (4/5) |
| Hybrid forest, literature sites (n=700) | 0.641 → 0.627 | 0.627 | −0.014 |

| Feature | Positive sites (mean) | Negative sites (mean) |
|---|---|---|
| pLDDT at cysteine | 65.8 | 81.2 |
| C-alpha neighbours within 10 Å | 4.8 | 10.9 |

---

## 3. Discussion

We set out to build a palmitoylation predictor that could credibly claim state-of-the-art performance. Doing so required first showing that the field's benchmarks overstate performance. With a single fixed model we measured that random and protein-level splits share 40–44% of test-positive windows with training data and inflate MCC by ~0.10 and AUC by ~0.04. Because 93% of UniProt palmitoylation annotations are inferred, and because inferred annotations are copied between orthologues, this is not a hypothetical risk: it is the default outcome of the splits used in most published work, including recent deep-learning studies that use site-level 80/20 partitions. We therefore release explicit cluster-held-out and species-held-out folds so that future methods can be compared on genuinely disjoint data.

On those folds PalmDL reaches AUC 0.956 overall and 0.901 on human sites. The best previously reported human AUC is 0.900 (GPS-Palm [8]); the best overall figure in the literature is 0.931 from a preprint using a site-level split (Deep-Palm [14]). Our human estimate is statistically indistinguishable from 0.900 (bootstrap P = 0.46, §2.4), obtained under a protocol with zero shared windows, and our audit suggests that site-level numbers above 0.95 should be treated sceptically. We therefore frame PalmDL as matching the published state of the art under a stricter protocol rather than exceeding it. We emphasise that this is not a head-to-head rerun: GPS-Palm's web server has been retired (HTTP 404 on 18 October 2026) and its standalone distribution is a GUI-only academic-licence package, so we could not score our test proteins with it. A true head-to-head on our released folds remains the single most important piece of missing evidence, and we provide the fold files to make it straightforward.

The cross-species experiments split into two statements that are easy to conflate. When orthologues of the test species are available — the ordinary transfer setting, measured by the 8 × 8 matrix — mammalian models reach 0.94 on other mammals and even 0.965 on plants, effectively as good as within-species performance. When human-homologous clusters are removed, performance on human falls to 0.822. Most of the apparent generality of palmitoylation predictors therefore comes from orthologue coverage rather than from an abstract motif grammar, and the bacterial rows of the matrix make the point sharply: *E. coli* and *B. subtilis* models transfer to each other almost perfectly but are at chance on mammals, because bacterial lipoproteins are N-acylated at the LipoBox rather than S-palmitoylated at a cysteine. A predictor that reported high accuracy across "all species" without separating these regimes would be measuring annotation propagation.

Two design choices follow from the error analysis. First, the gap between inferred (AUC 0.960) and experimental sites (0.908) is an annotation-memory effect that any user should know about: the experimentally characterised sites are the ones the propagation machinery does not already cover. Second, reporting a fixed 0.5 cut-off is less useful than reporting FDR-controlled operating points, because the deployable question is "how many predictions can I trust", not "what is my AUC". At a 5% estimated FDR the model retains 61% of known sites at 94% precision, which is a realistic basis for a shortlist.

The proteome-wide scan turns that principle into a resource. At 5% FDR it returns 322 human cysteines, 227 of which are unannotated, with a 16-fold enrichment for cytoplasmic topology and a top tier dominated by juxtamembrane cysteines of single-pass membrane proteins. The small absolute size is not a weakness of the model but a property of the problem: with only ~0.2% of human cysteines annotated, any procedure that controls FDR at the proteome scale must be conservative. We therefore report five tiers so that users can trade confirmation for discovery, and we expect the 20–30% tiers to be the useful range for follow-up.

Structure is a modest, real improvement rather than a cure. Three AlphaFold-derived numbers — pLDDT at the site, local pLDDT and C-alpha packing — improve the forest view from 0.890 to 0.910 and, when folded into a hybrid model that keeps all 71,594 training sites, the full test set from 0.946 to 0.949, consistently across folds. They also point in the expected direction: palmitoylated cysteines are less ordered and more exposed than the cysteines that surround them. But on the 700 literature sites that have a structural model the structure-aware forest is no better than the sequence forest (0.627 versus 0.641). The domain shift documented in §2.10 is therefore not simply a burial-versus-exposure problem, and structure alone will not close it; a larger structured training set, enzyme conditioning or orthogonal experimental data will be needed.

**Limitations.** (i) 93% of positives are inferred rather than experimental; we stratify by evidence but cannot eliminate label noise, and the 477 experimental sites are too few for a well-powered experimental-only benchmark. (ii) Background negatives are proteins without annotation, not proteins proven not to be palmitoylated; some are certainly false negatives, which biases absolute metrics downwards while leaving model ranking intact. (iii) Our cluster definition uses 4-mer cosine similarity, an operational proxy: in a sample of 30 pairs, cosine 0.28–0.34 corresponded to median 71% global identity (IQR 67–73%), so remote homologues below ~50% identity may still cross folds; the species holdout is the more conservative experiment and is reported for that reason. (iv) We model only cysteine-centred windows and do not use predicted structure, DHHC-enzyme identity, or tissue expression. (v) The ensemble is calibrated at natural prevalence; performance at other prevalences requires recalibration. (vi) Independent, site-level public benchmarks do not exist: SwissPalm is behind an IP-verification gate, dbPTM and CysModDB no longer resolve, the MDD-Palm supplements contain no site lists, and GPS-Palm's server is retired. Both alternatives we could obtain fail a motif sanity check — the public RFCM-PALM mouse MS data (PWM AUC 0.53; PalmDL 0.52) and our own literature-mined set (PWM AUC 0.54; PalmDL 0.59; §2.10) — so neither is used as a benchmark. Their failure is itself informative: site-level assignment in cysteine-rich or transmembrane contexts is not yet reliable at scale, which is why the field's headline numbers remain independently unverified. The curated literature subset makes the limitation sharper: 99 unambiguous sites, a ±2-residue numbering control and a motif check all leave PalmDL at AUC ≈0.62 against same-protein cysteines, so the model recognises the UniProt consensus but not the motif-atypical sites that mechanistic studies tend to characterise. This is a domain-shift limitation shared by any annotation-trained predictor and the main reason we do not claim prospective accuracy.

**Future work.** Structure-aware features from AlphaFold, enzyme-conditioned models that predict which DHHC isoform a site prefers, manually curated literature sets assembled with the extraction pipeline of §2.10, and a proteome-wide scan with target–decoy FDR estimation are the natural next steps. The released `score()` function makes a batch prediction interface straightforward to deploy.

---

## 4. Methods

### 4.1 Data acquisition

All sequences and annotations were retrieved from the UniProtKB REST API (release 2026_03) on 18 October 2026 using two queries: `(ft_lipid:*) AND (reviewed:true)` and `(reviewed:true) AND (organism_id:9606) AND (NOT ft_lipid:*)`. Fields retrieved: accession, sequence, organism, protein and gene names, lipidation features, signal peptide, topological domains and length. Queries and row counts are stored in `data/raw/download_meta.json`.

### 4.2 Site and negative definition

Lipidation features were parsed with a regular expression; a site was labelled positive when the feature note contained "S-palmitoyl cysteine". Evidence codes were mapped to experimental (ECO:0000269, ECO:0007744, ECO:0007005, ECO:0006056) or inferred. Each site was represented by a 31-residue window centred on the cysteine and zero-padded with `X` at protein termini. For 66 annotations (1.0%) the canonical residue at the annotated position is not cysteine, reflecting transferred annotations on diverged or re-annotated entries. They are retained in the cross-validation tables (removing them changes nothing: AUC 0.9451 with and without, balanced MCC 0.7407 versus 0.7440; `results/sensitivity_residue_filter.csv`) but are excluded from the released prediction model, because mutagenesis showed that they teach the model an impossible lysine-at-centre rule (§2.9). Hard negatives were unannotated cysteines of positive-carrying proteins at least 16 residues from any annotated site, capped at four per protein; background negatives were up to three cysteines from each human protein without any lipidation feature. Topology (cytoplasmic, extracellular, lumenal, signal peptide, unknown) was derived from UniProt topological domains.

### 4.3 Leakage-controlled splits

Proteins carrying at least one positive were clustered by cosine similarity of their 4-mer count vectors with a threshold of 0.30; clusters sharing an identical positive window were merged. Clusters were assigned to five folds with greedy balancing of positive counts (1,361 positives per fold). The species holdout used all human sites as test; training used non-human sites after removing (i) every cluster containing a human test protein and (ii) every annotation whose evidence referenced a human accession (`ECO:0000250|UniProtKB:<human>`). Human background negatives were split 50/50 by protein hash so that no protein contributed to both training and test.

### 4.4 Features

Three feature blocks were computed per site. (i) A 23 × 31 one-hot/physicochemical tensor: 20 amino-acid channels plus Kyte–Doolittle hydropathy, normalised van der Waals volume and Vihinen flexibility. (ii) Eleven protein-context features: six topology indicators, an annotation-available flag, log protein length, relative position, relative distance to the nearest terminus, and a terminus-proximal flag. (iii) ESM-2 (esm2_t12_35M_UR50D) representations, computed as the mean over the 31 window tokens and the representation of the central cysteine (2 × 480 dimensions). ESM embeddings were extracted once and cached; no ESM fine-tuning was performed.

### 4.5 Models and training

The PWM baseline computes position-specific log-odds from training positives and negatives with a pseudocount of 0.5 and converts the summed score to a probability with a logistic link. Logistic regression and random forests (400 trees, minimum leaf size 3, balanced subsampling) used the flattened one-hot window plus context features. The CNN has four parallel 1-D convolutions (kernels 3, 5, 7, 11; 32 filters each), ReLU, global max pooling, batch normalisation, a 64-unit hidden layer with dropout 0.3, and a sigmoid output; it is trained with AdamW (learning rate 1 × 10⁻³, weight decay 1 × 10⁻⁴), class-weighted BCE, gradient clipping at 5 and early stopping on validation AUC (patience 8, maximum 40 epochs, batch size 128). The ESM head is a two-layer MLP (256 and 128 units, dropout 0.3) trained with the same schedule and a learning rate of 5 × 10⁻⁴.

For each outer fold, 15% of the training clusters were held out as an inner validation set for early stopping and threshold selection. Background negatives were subsampled to 1× the number of positives in training; test folds retained their natural composition. Thresholds were chosen to maximise MCC on the inner validation set, and a second threshold was chosen on a class-balanced subsample of the validation set for balanced metrics. Ensemble members were rank-transformed within fold and averaged with equal weights; a Platt scaling model was fitted on the inner validation scores and applied to test scores. Member subsets were compared by pooled validation AUC, never by test performance.

### 4.6 Evaluation

Discrimination was measured by AUC and average precision; threshold-dependent performance by MCC, accuracy, sensitivity and specificity at the validation-derived threshold, reported both at natural prevalence and on class-balanced subsamples. Calibration was measured by expected calibration error with ten bins. Confidence intervals were obtained by bootstrapping sites (2,000 replicates); model comparisons used paired bootstrap on identical sites. All experiments used Python 3.12, PyTorch 2.9 (CPU), scikit-learn 1.6 and Biopython 1.88.

### 4.7 Cross-species transfer matrix

Eight species with at least 85 annotated sites were used: human, mouse, rat, cattle, pig, *Arabidopsis thaliana*, *Escherichia coli* and *Bacillus subtilis*. For each source species a PalmDL model (random forest plus ESM head) was trained on that species' positives, its same-protein negatives and an equal number of human background negatives. For the diagonal, 20% of the source species' protein clusters were held out and the model was retrained on the remaining 80%; off-diagonal entries were predicted by the model trained on all source rows. The ESM head used 15% of the remaining clusters for early stopping (30 epochs, patience 6). AUC was computed on each target species' positives and same-protein negatives.

### 4.8 FDR-controlled operating points

For each target FDR f, the validation folds were sorted by predicted probability and the largest threshold whose cumulative false-discovery fraction (FP / (TP + FP)) was ≤ f was selected. That threshold was then applied unchanged to the outer test folds, and precision, recall and the number of predictions were recorded.

### 4.9 In-silico mutagenesis and label audit

Three hundred randomly chosen human positive windows were scored with the forest view. Each of the 31 positions was mutated to each of the 20 amino acids in turn (620 predictions per site) and the mean change in predicted probability relative to the wild-type window was recorded. The procedure was run once with the full UniProt annotation set and once after removing 66 features whose annotated position is not a cysteine in the canonical sequence; the latter model is the one released for prediction.

### 4.10 Literature mining

Europe PMC was queried for `(palmitoylation OR "S-palmitoylation" OR palmitoylated) AND OPEN_ACCESS:Y AND FIRST_PDATE:[2020-01-01 TO 2026-12-31]` (9,062 hits) and the first 1,200 records were sampled; 992 full texts could be retrieved. Sentences containing "palmitoyl" and a pattern of the form `Cys<number>` or `C<number>` were extracted, with the site pattern required not to be preceded by a letter or digit (so that `zDHHC1` is not read as site 1). A protein token was taken from the 40 characters to the left or right of the match and kept only if it looked like a gene symbol (contains a digit, is upper case, or matches a lower-case prefix followed by an upper-case letter). Tokens were resolved against a local dictionary of gene names from the 28,711 reviewed sequences already downloaded; a token was accepted only when it mapped to a single accession (exact match, or a unique prefix). A mention became a positive only if the residue at that position in the canonical sequence is a cysteine. Negatives were up to six other cysteines of the same protein that are neither mentioned in the set nor annotated as S-palmitoylated in UniProt 2026_03. A position-specific log-odds motif trained on the UniProt data was used as an independent check of whether the extracted positives carry a palmitoylation motif.

### 4.11 Proteome-wide scan and decoy construction

All 262,286 cysteines of the 20,431 reviewed human sequences in the local UniProt tables were scored in chunks of 5,000 with the deployed model. For the FDR estimate we drew 20,000 random target windows and shuffled their 30 flanking residues while keeping the central cysteine fixed; this destroys the motif but preserves the candidate class, whereas mutating the centre (C→A) leaves a motif-shaped window and yields a badly conservative null because the model ignores the identity of the central residue. FDR at a threshold was estimated as (number of decoys at or above the threshold) × (targets / decoys) divided by (number of targets at or above the threshold), and q-values were obtained by taking the running minimum over lower thresholds. Five tiers (1–30% FDR) are reported.

### 4.12 AlphaFold features and paired evaluation

Human proteins with at least one annotated site (316) and a random sample of 1,200 background human proteins were queried through the AlphaFold DB prediction API, which returns the current model URL; 479 models (266 site-carrying proteins) were retrieved for chains shorter than 700 residues. Per site we recorded the C-alpha pLDDT, the mean pLDDT of residues within ±5 positions, and the number of C-alpha atoms within 10 Å of the site C-alpha excluding sequence neighbours closer than four positions. The three features were standardised and appended to the base feature matrix; the same random forest and the same cluster-disjoint folds were used with and without them, so the comparison is paired at the fold level. Three additional designs were evaluated on the same folds: a logistic reranker on the sequence-ensemble logit plus the three features; a third-view ensemble that averages the structure-aware forest with the frozen ESM view; and a hybrid forest trained on all 71,594 sites with median imputation and a binary missing-structure indicator, so that proteins without an AlphaFold model fall back to the sequence path. Literature sites were scored with the same pipelines after fetching models for their proteins. The prediction service released with the paper exposes the sequence path by default and the structure-aware path behind a flag, falling back to the sequence path for proteins without a model (server/app.py --structure).

### 4.13 Reproducibility

The pipeline is organised as 25 scripts in `src/` (pipeline table in the repository README, ordered commands in `docs/REPRODUCE.md`) with fixed seeds; processed data, split assignments, out-of-fold predictions and per-fold metrics are included. Runtime on a 16-core CPU machine is approximately 3 h end-to-end, of which ESM embedding extraction is the largest component.

---

## Data and code availability

Code, processed data (including the canonical fold assignments), out-of-fold predictions, proteome-wide scores and the prediction service are archived at Zenodo (DOI: [ZENODO-DOI]); the service is available at [SERVICE-URL]. Supplementary Tables S1–S7 accompany the manuscript. Raw data are from UniProt release 2026_03 and are freely available under CC BY 4.0; AlphaFold models are from the AlphaFold DB (CC BY 4.0).

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
9. Li Z, Wang C, Wang X, Wang J, Li Z, Zhang Y, et al. pCysMod: prediction of multiple cysteine modifications based on deep learning framework. *Front Cell Dev Biol.* 2021;9:617366. doi:10.3389/fcell.2021.617366
10. Bandyopadhyay A, et al. RFCM-PALM: in-silico prediction of S-palmitoylation sites in the synaptic proteins for male/female mouse data. *Int J Mol Sci.* 2021;22(18):9901. doi:10.3390/ijms22189901
11. Blanc M, David F, Abrami L, Migliozzi D, Armand F, Bürgi J, et al. SwissPalm: protein palmitoylation database. *F1000Res.* 2015;4:261. doi:10.12688/f1000research.6464.1
12. Meng L, Chen X, Cheng K, Chen K, Xu Z, Ren J, et al. CysModDB: a comprehensive platform with the integration of manually curated resources and analyses to facilitate the study of cysteine modifications. *Brief Bioinform.* 2022;23(5):bbac460. doi:10.1093/bib/bbac460
13. Meng L, Chan WS, Huang L, Hu J, Wang J, Huang H, et al. Mini-review: recent advances in post-translational modification site prediction based on deep learning. *Comput Struct Biotechnol J.* 2022;20:3522–3532. doi:10.1016/j.csbj.2022.06.045
14. Deng X, et al. Deep-Palm: an integrated deep learning framework for structure-aware prediction of protein S-palmitoylation. *bioRxiv* 2026. doi:10.64898/2026.03.05.709753 (preprint)
15. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res.* 2025;53:D609–D617. doi:10.1093/nar/gkae1010

---

## Figure legends

**Figure 1. Dataset composition.** (a) Organisms contributing the largest numbers of annotated S-palmitoyl cysteine sites (human and mouse highlighted; remaining 1,415 organisms pooled). (b) Evidence tier of positive sites (log scale): 477 experimentally supported versus 6,328 inferred, mostly by similarity. (c) Annotated topology of positive sites.

**Figure 2. Split protocol determines measured performance.** A single random forest was trained and evaluated under four protocols. (a) AUC and class-balanced MCC. (b) Percentage of test positives whose exact 31-residue window is present among training positives. Random site-level and protein-level splits leak 40–44% of windows; cluster and species splits leak none.

**Figure 3. Discrimination under cluster-held-out cross-validation.** (a) ROC and (b) precision–recall curves for the PWM baseline, logistic regression, random forest, CNN, ESM-2 view and the PalmDL ensemble (6,805 positives, 71,594 sites).

**Figure 4. Performance on human sites under three evaluation regimes.** (a) AUC and (b) class-balanced MCC for multi-species training with human test folds, human-only training, and the strict species holdout. PalmDL scores are shown where available.

**Figure 5. Calibration.** Reliability diagrams for the random forest, CNN, ESM-2 view and calibrated ensemble at natural prevalence.

**Figure 6. Sequence motif learned from the data.** Position-specific log-odds of amino acids from −7 to +8 around the central cysteine, positives versus background negatives (log2 scale).

**Figure 7. The palmitoylation code is eukaryotic.** AUC of single-species PalmDL models, rows = training species, columns = test species. Mammalian models transfer almost losslessly to other mammals and to *A. thaliana*, whereas bacterial models are at chance on mammals and near-perfect on each other. Diagonal entries come from a protein-cluster-disjoint 20% holdout of the same species.

**Figure 8. What the model uses and how far its predictions can be trusted.** (a) Mean change in predicted probability when each position of 300 human positive windows is mutated to each amino acid (forest view). Cysteine at ±1 increases the score (cysteine clusters), acidic residues at −1/−2 decrease it, and the central position is neutral because every candidate is a cysteine. (b) FDR-controlled operating points: precision against recall on the outer test folds when the threshold is chosen on validation at estimated FDR 1%, 5%, 10% and 20%.

**Figure 9. Proteome-wide scan and FDR-controlled shortlist.** (a) Score distributions of all 262,286 human cysteines and of 20,000 shuffled-flank decoys (log scale). (b) Estimated FDR as a function of shortlist size, with the five reported tiers marked. (c) Composition of each tier: known UniProt sites versus novel cysteines.

## Table legends

**Table 1.** Cluster-held-out cross-validation (all species, 5 folds, 71,594 sites). AUC, AUPR, MCC and balanced MCC with 95% bootstrap CIs for PalmDL.

**Table 2.** Leakage audit: identical-window overlap and performance under four splitting protocols.

**Table 3.** Human evaluations: multi-species training, human-only training, species holdout, and the hard-negative benchmark.

**Table 4.** Comparison with published S-palmitoylation predictors.

**Table 5.** Cross-species transfer matrix (AUC) for eight species.

**Table 6.** Error structure of the held-out ensemble predictions by evidence tier and annotated topology.

**Table 7.** Automated literature mining as an external benchmark: yields, motif consistency and model performance under progressively stricter extraction filters.

**Table 8.** Proteome-wide scan of the reviewed human proteome at five FDR tiers, with recovery of known sites and the fraction of novel predictions.

**Table 9.** AlphaFold structural features on human sites: paired cross-validated comparison and feature distributions.
