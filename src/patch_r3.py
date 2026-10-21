"""Third-round revisions: Table 6/ECE from nested, matrix text, extrapolation limits."""
from __future__ import annotations

from pathlib import Path

MD = Path(__file__).resolve().parents[1] / "manuscript" / "manuscript.md"
t = MD.read_text(encoding="utf-8")


def rep(old: str, new: str) -> None:
    global t
    if old not in t:
        print("MISS:", old[:70])
        return
    t = t.replace(old, new)


# --- Table 6 (pooled nested) ---
rep("| Experimentally supported sites | 455 | 0.907 | 0.176 |", "| Experimentally supported sites | 455 | 0.902 | 0.176 |")
rep("| Inferred sites | 3,969 | 0.947 | 0.682 |", "| Inferred sites | 3,969 | 0.941 | 0.664 |")
rep("| Cytoplasmic | 2,445 | 0.887 | 0.829 |", "| Cytoplasmic | 2,445 | 0.858 | 0.808 |")
rep("| Extracellular | 27 | 0.855 | 0.375 |", "| Extracellular | 27 | 0.828 | 0.361 |")
rep("| Topology unknown | 1,801 | 0.898 | 0.501 |", "| Topology unknown | 1,801 | 0.900 | 0.512 |")

# --- §2.8 text ---
rep("Discrimination remains better for inferred than for experimentally supported sites (AUC 0.947 versus 0.907).",
    "Discrimination remains better for inferred than for experimentally supported sites (pooled AUC 0.941 versus 0.902; each tier uses its positives plus all 62,842 negatives).")
rep("Performance is stable across annotated topologies (cytoplasmic 0.887, extracellular 0.855, topology unknown 0.898), with the lumenal class too small to interpret (14 positives).",
    "Pooled performance is stable across annotated topologies (cytoplasmic 0.858, extracellular 0.828, topology unknown 0.900); the extracellular estimate draws on three of five folds with both classes and the lumenal class (14 positives) is not interpretable. These pooled values differ from the unweighted fold means in Table 1, which are the primary estimates.")

# --- ECE ---
rep("Probabilities are Platt-calibrated (mean fold ECE 0.025).",
    "Probabilities are Platt-calibrated (mean fold ECE 0.024; pooled ECE 0.018).")
rep("Probabilities from the calibrated ensemble are reliable at natural prevalence (mean fold ECE 0.025, pooled 10-bin ECE 0.021; Figure 5)",
    "Probabilities from the calibrated ensemble are reliable at natural prevalence (mean fold ECE 0.024, pooled 10-bin ECE 0.018; Figure 5)")

# --- counts 449/3975 -> 455/3969 and 89.9 -> 89.7 ---
rep("449 experimentally supported; 786 organisms", "455 experimentally supported; 786 organisms")
rep("455 experimentally supported; 786 organisms; 557 human sites", "455 experimentally supported; 786 organisms; 557 human sites")
rep("449 supported by experimental evidence codes and 3,975 inferred (mostly by similarity, 89.9%)",
    "455 supported by experimental evidence codes and 3,969 inferred (mostly by similarity, 89.7%)")
rep("449 experimental", "455 experimental")
rep("3,975 inferred", "3,969 inferred")

# --- SOTA / cross-paper claims ---
rep("We show that PalmDL exceeds the best previously reported human AUC on a benchmark with zero shared windows, and we release the split definitions",
    "We report PalmDL's human AUC on a benchmark with zero shared windows alongside published values obtained on different data, and we release the split definitions")
rep("Our human estimate is marginally above 0.900", "Our human estimate is 0.910")
rep("PalmDL's 0.938 is obtained with zero shared windows under nested member selection, and its human estimate, 0.910, is reported for context rather than as a paired superiority claim.",
    "PalmDL's 0.938 is obtained with zero shared windows under nested member selection, and its human estimate, 0.910, is reported for context rather than as a paired superiority claim, because the published values come from different data, negatives and splits.")

# --- FDR naming ---
rep("### 2.8 Error structure and FDR-controlled shortlists", "### 2.8 Error structure and validation-target operating points")
rep("and calibrated scores support FDR-controlled shortlists", "and calibrated scores support validation-target candidate shortlists")
rep("For practical use we derived, **within each fold**, the lowest threshold whose complete score groups satisfy an empirical validation false-discovery target",
    "For practical use we derived, **within each fold**, the lowest threshold whose complete score groups satisfy an empirical validation-target operating point")
rep("### 2.11 A proteome-wide, FDR-controlled candidate list", "### 2.11 A proteome-wide candidate list with empirical q-values")
rep("we applied the forest view of the released model to every cysteine of the reviewed human proteome — **262,286 sites in 19,814 proteins** — and controlled the false-discovery rate with 20,000 shuffled-flank decoys",
    "we applied the forest view of the released model to every cysteine of the reviewed human proteome — **262,286 sites in 19,814 proteins** — and estimated empirical q-values against 20,000 shuffled-flank decoys")
rep("**Table 8. Proteome-wide scan of the reviewed human proteome (262,286 cysteines) at five FDR tiers.**",
    "**Table 8. Proteome-wide scan of the reviewed human proteome (262,286 cysteines) at five empirical candidate tiers.**")

# --- literature / generalization ---
rep("The curated subset behaves like the full set, so extraction noise is not the explanation.",
    "The curated subset also has little overlap with the candidate tiers (0/1/2 at 5%/20%/30%), consistent with the full set; extraction noise alone is not a sufficient explanation, but the curated set has no independent AUC of its own.")
rep("The honest conclusion is a domain-shift limitation that applies to any annotation-trained predictor, including PalmDL:",
    "For PalmDL, the honest conclusion is a domain-shift limitation:")
rep("It quantifies the statement in the Introduction: the field's headline numbers have never been checked on data that are simultaneously independent, site-level and motif-consistent, and automated literature mining does not yet close that gap.",
    "It quantifies the statement in the Introduction for our own model: its headline numbers have not been checked on data that are simultaneously independent, site-level and motif-consistent, and automated literature mining does not yet close that gap. A low motif score does not establish that a literature site is wrong; it shows that the site lies outside the motif distribution the model learned.")

# --- structure limits ---
rep("The domain shift documented in §2.10 is not simply a burial-versus-exposure problem, and structure alone will not close it.",
    "No structural model was available for the literature sites in this revision, so whether structure would help there remains untested rather than refuted.")

# --- matrix text corrections ---
rep("Most mammalian pairs fall between 0.83 and 0.95 (mean 0.896), with Pig→Human (0.834) the clear exception, 0.11 below the human diagonal and 0.14 below the pig diagonal.",
    "Most mammalian pairs fall between 0.83 and 0.95; the arithmetic mean of the 18 observable mammalian off-diagonal cells is 0.901. Pig→Human (0.834) is the clear exception, 0.084 below the human diagonal. Pig's own diagonal is empty (fewer than 20 eligible positives), so no diagonal comparison is made for that row.")
rep("The plant model reaches 0.74–0.78 on mammals, whereas mammalian models reach 0.96–0.98 on *A. thaliana*.",
    "The plant model reaches 0.739–0.785 on mammals, whereas mammalian models reach 0.96–0.98 on *A. thaliana*.")
rep("Two consequences follow. First, claims of a single shared eukaryotic code are not supported by this design: the pig exception and the plant asymmetry are as informative as the mammalian mean.",
    "Two consequences follow. First, claims of a single shared eukaryotic code are not supported by this design: the pig exception and the plant asymmetry are as informative as the mammalian mean. For the human source only, the early-stopping split also contained 6,298 background negatives (the other sources use same-protein negatives), so the human row's model-selection condition is not fully matched; this is disclosed rather than corrected in this revision.")
rep("With orthologues available, most mammalian pairs transfer at 0.83–0.95, but Pig→Human is 0.834 and the plant→mammal direction is 0.74–0.78;",
    "With orthologues available, most mammalian pairs transfer at 0.83–0.95 (18-cell mean 0.901), but Pig→Human is 0.834 and the plant→mammal direction is 0.739–0.785;")
rep("its human estimate, 0.910, is reported for context", "its human estimate, 0.910, is reported for context")

# Table 5 note
rep("One model per source species is trained on its positives and same-protein negatives (no human background); 20% of source clusters are held out for the diagonal and another 15% for early stopping.",
    "One model per source species is trained on its positives and same-protein negatives (no human background); 20% of source clusters are held out for the diagonal and another 15% for early stopping. For the human source the early-stopping split also contained 6,298 background negatives, unlike the other rows; the 18-cell mammalian off-diagonal mean is 0.901.")

# Figure 1 legend counts and percentage
rep("(b) Evidence tier of positive sites (log scale): 449 experimentally supported versus 3,975 inferred, mostly by similarity.",
    "(b) Evidence tier of positive sites (log scale): 455 experimentally supported versus 3,969 inferred, mostly by similarity (89.7%).")

# Figure 8 panel text
rep("(b) FDR-controlled operating points: precision against recall on the outer test folds when the threshold is chosen on validation at estimated FDR 1%, 5%, 10% and 20%.",
    "(b) Validation-target operating points: precision against recall on the outer test folds when the threshold is chosen on validation at empirical targets of 1%, 5%, 10% and 20%; observed held-out false-discovery proportions are 4.7%, 10.0%, 15.1% and 24.1%.")

# Data availability: no "or later"
rep("available at https://github.com/opus940810/palm-dl (commit 758952f or later).",
    "available at https://github.com/opus940810/palm-dl (release tag v2.1; the exact commit, file hashes and reproduction commands are listed in docs/REPRODUCE.md).")

# References 10 and 12
rep("10. Bandyopadhyay A, et al. RFCM-PALM:", "10. Bandyopadhyay SS, et al. RFCM-PALM:")
rep("12. Meng Y, Zhang L, Zhang L, Wang Z, Wang X, Li C. CysModDB: a comprehensive platform with the integration of manually curated resources and analyses to facilitate the study of cysteine modifications.",
    "12. Meng Y, Zhang L, Zhang L, Wang Z, Wang X, Li C. CysModDB: a comprehensive platform with the integration of manually curated resources and analysis tools to facilitate the study of cysteine modifications.")

MD.write_text(t, encoding="utf-8")
print("third-round manuscript patch applied")
