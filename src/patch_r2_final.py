"""Final text fixes: references, availability, headings, aggregation, Table 7."""
from __future__ import annotations

from pathlib import Path

MD = Path(__file__).resolve().parents[1] / "manuscript" / "manuscript.md"
t = MD.read_text(encoding="utf-8")


def rep(old: str, new: str) -> None:
    global t
    if old not in t:
        print("MISS:", old[:60])
        return
    t = t.replace(old, new)


rep("### 2.7 The palmitoylation code is shared across mammals and partly with plants",
    "### 2.7 Cross-species transfer is heterogeneous")

rep("Training on human data alone is selected as RF+ESM in four folds and ESM alone in one (mean 0.904).",
    "Training on human data alone is selected as RF+ESM in four folds and CNN+ESM in one (mean 0.901).")

rep("9. Li Z, Wang C, Wang X, Wang J, Li Z, Zhang Y, et al. pCysMod:",
    "9. Li S, Yu K, Wu G, Zhang Q, Wang P, Zheng J. pCysMod:")
rep("12. Meng L, Chen X, Cheng K, Chen K, Xu Z, Ren J, et al. CysModDB: a comprehensive platform with the integration of manually curated resources and analyses to facilitate the study of cysteine modifications. *Brief Bioinform.* 2022;23(5):bbac460.",
    "12. Meng Y, Zhang L, Zhang L, Wang Z, Wang X, Li C. CysModDB: a comprehensive platform with the integration of manually curated resources and analyses to facilitate the study of cysteine modifications. *Brief Bioinform.* 2022;23(6):bbac460.")
rep("13. Meng L, Chan WS, Huang L, Hu J, Wang J, Huang H, et al. Mini-review:",
    "13. Meng L, Chan W, Huang L, Liu L, Chen X, Zhang W. Mini-review:")
rep("14. Deng X, et al. Deep-Palm:", "14. Deng M, et al. Deep-Palm:")

rep("Code, processed data (including the canonical fold assignments), out-of-fold predictions, proteome-wide scores and the prediction service are archived at Zenodo (DOI: 10.5281/zenodo.23455873); the service is available at [SERVICE-URL]. Supplementary Tables S1–S7 accompany the manuscript.",
    "The corrected v2 code, processed data (including the canonical fold assignments), out-of-fold predictions and proteome-wide scores are available at https://github.com/opus940810/palm-dl (commit 758952f or later). The earlier v1.0.0 release is archived at Zenodo (DOI 10.5281/zenodo.23455873); it predates the label correction and must not be used for the results reported here, and a corrected Zenodo version will be deposited before publication. The prediction service is available at [SERVICE-URL]. Supplementary Tables S1–S7 accompany the manuscript.")

rep("| All verified sites | 133 | 0.537 | 0.588 |",
    "| All verified sites | 133 | 0.537 | 0.585† |")
rep("Comparisons are against other cysteines of the same proteins. PWM AUC is the enrichment of a position-specific log-odds motif trained on the UniProt data; if the extracted positives were faithful, this value should approach the model AUC. The manually curated subset (99 unambiguous mentions in 82 proteins) is used for the overlap analysis rather than for a separate filter row.",
    "†Weighted within-protein AUC over 104 protein-level tests; the pooled novel-site AUC against same-protein cysteines is 0.578 and the filter-table value for the same 133 sites is 0.588. Comparisons are against other cysteines of the same proteins. The manually curated subset (99 unambiguous mentions in 82 proteins) is used for the overlap analysis, not as a separate filter row.")

rep("Confidence intervals were obtained by bootstrapping sites (2,000 replicates); model comparisons used paired bootstrap on identical sites.",
    "Confidence intervals were obtained by bootstrapping sites (2,000 replicates) and, because sites within a protein cluster are not independent, by resampling protein clusters within folds (1,000 replicates); the wider cluster interval is the one quoted in the abstract. Model comparisons used paired cluster bootstrap on identical rows.")

MD.write_text(t, encoding="utf-8")
print("final text fixes applied")
