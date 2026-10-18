"""Build the supplementary document (Tables S1-S7)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from docx import Document
from docx.shared import Pt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "manuscript" / "PalmDL_supplementary.docx"
FIG = ROOT / "figures"


def add_table(doc, title: str, header: list[str], rows: list[list], note: str = "", sizes=(7.5, 7.0)):
    doc.add_heading(title, level=2)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        r = c.paragraphs[0].add_run(h)
        r.bold = True
        r.font.size = Pt(sizes[1])
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row[: len(header)]):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(str(v))
            run.font.size = Pt(sizes[1])
    if note:
        p = doc.add_paragraph()
        r = p.add_run(note)
        r.italic = True
        r.font.size = Pt(7.5)


def protein_names() -> dict:
    names = {}
    for f in ("uniprot_negatives.tsv", "uniprot_lipid_reviewed.tsv"):
        d = pd.read_csv(ROOT / "data" / "raw" / f, sep="\t", dtype=str).fillna("")
        names.update(dict(zip(d["Entry"], d["Protein names"])))
    return names


def main() -> None:
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(9)
    doc.add_heading("Supplementary Material", level=1)
    p = doc.add_paragraph()
    r = p.add_run("PalmDL: a leakage-aware deep learning framework and benchmark for protein "
                  "S-palmitoylation site prediction")
    r.bold = True
    r.font.size = Pt(11)
    doc.add_paragraph("This supplement contains the sensitivity analyses, per-fold metrics, full "
                      "leakage table, curated literature set, proteome-scan tiers, structure "
                      "comparisons and the statistical comparison against published values.")

    # S1 sensitivity
    s1 = pd.read_csv(RES / "sensitivity_residue_filter.csv")
    add_table(doc, "Table S1. Sensitivity to 66 annotations whose canonical residue is not cysteine.",
              ["Dataset", "n sites", "n positives", "AUC", "AUC SD", "Balanced MCC"],
              [[r.label, int(r.n), int(r.n_pos), f"{r.auc:.4f}", f"{r.auc_sd:.4f}", f"{r.bal_mcc:.4f}"]
               for r in s1.itertuples()],
              "The same folds and model; only the 66 non-cysteine annotations are removed.")

    # S2 per-fold CV metrics
    frames = []
    for f, tag in (("metrics_cv_base.csv", None), ("metrics_cv_cnn.csv", None), ("metrics_cv_esm.csv", None),
                   ("metrics_cv_ensemble.csv", None)):
        p = RES / f
        if p.exists():
            d = pd.read_csv(p)
            d = d[d.protocol == "test_natural"]
            frames.append(d)
    if frames:
        m = pd.concat(frames)
        agg = (m.groupby("model")[["mcc", "auc", "aupr", "acc", "sn", "sp"]]
               .agg(["mean", "std"]).round(4))
        rows = []
        for model, r in agg.iterrows():
            rows.append([model] + [f"{r[(c, 'mean')]:.3f} ± {r[(c, 'std')]:.3f}"
                                   for c in ["auc", "aupr", "mcc", "acc", "sn", "sp"]])
        add_table(doc, "Table S2. Per-fold cross-validation metrics (mean ± SD over five folds).",
                  ["Model", "AUC", "AUPR", "MCC", "Accuracy", "Sn", "Sp"], rows,
                  "Cluster-held-out folds, all species, natural prevalence; thresholds chosen on "
                  "inner validation folds.")

    # S3 leakage ladder
    s3 = pd.read_csv(RES / "leakage_ladder.csv")
    add_table(doc, "Table S3. Split-protocol audit (same random forest and features).",
              ["Protocol", "Train pos.", "Test sites", "Test pos.", "AUC", "MCC", "Balanced MCC",
               "Shared windows"],
              [[r.protocol, int(r.n_train_pos), int(r.n_test), int(r.n_test_pos), f"{r.auc:.3f}",
                f"{r.mcc_nat:.3f}", f"{r.mcc_bal:.3f}", f"{r.leaked_pos_windows * 100:.1f}%"]
               for r in s3.itertuples()],
              "'Shared windows' is the fraction of test positives whose exact 31-residue window "
              "occurs among training positives.")

    # S4 curated literature sites
    cur = pd.read_csv(RES / "literature" / "curated_sites.csv")
    cur = cur[cur.curated]
    names = protein_names()
    rows = [[i, r.acc, f"C{r.pos}", r.token, str(names.get(r.acc, ""))[:60], r.pmcid]
            for i, r in enumerate(cur.itertuples(), 1)]
    add_table(doc, "Table S4. Curated literature sites (99 of 136 extracted mentions).",
              ["#", "Accession", "Site", "Token", "Protein", "PMCID"], rows,
              "Accepted when the sentence attributes palmitoylation to that protein at that "
              "cysteine, or a cysteine mutant loses palmitoylation; rejected otherwise (§2.10).",
              sizes=(6.5, 6.0))

    # S5 scan tiers + enrichment
    fdr = pd.read_csv(RES / "scan" / "fdr_summary.csv")
    add_table(doc, "Table S5a. Proteome-wide scan at five FDR tiers.",
              ["FDR", "Sites", "Proteins", "Known recovered", "Recall", "Novel", "Threshold"],
              [[f"{r.fdr:.0%}", int(r.n_predictions), int(r.n_proteins), int(r.known_recovered),
                f"{r.known_recall:.1%}", int(r.novel), f"{r.threshold:.3f}"] for r in fdr.itertuples()],
              "262,286 cysteines in 19,814 reviewed human proteins; shuffled-flank decoys.")
    ce = pd.read_csv(RES / "scan" / "compartment_enrichment.csv")
    add_table(doc, "Table S5b. Topology enrichment of the 5% FDR shortlist.",
              ["Compartment", "Share in shortlist", "Share in all cysteines", "Enrichment"],
              [[r.compartment, f"{r.share_shortlist:.1%}", f"{r.share_all:.1%}", f"{r.enrichment:.1f}×"]
               for r in ce.itertuples()],
              "Topology from UniProt; 'unknown' dominates the background and is depleted.")

    # S6 structure comparisons
    h = pd.read_csv(RES / "structure" / "hybrid_cv.csv")
    rows = []
    for fold, g in h.groupby("fold"):
        b = g[g.model == "base"].iloc[0]
        s = g[g.model == "hybrid"].iloc[0]
        rows.append([int(fold), f"{b.auc_all:.3f}", f"{s.auc_all:.3f}",
                     f"{b.auc_structured:.3f}", f"{s.auc_structured:.3f}"])
    add_table(doc, "Table S6. Hybrid structure-aware forest, paired per fold.",
              ["Fold", "Base AUC (all)", "Hybrid AUC (all)", "Base AUC (structured)",
               "Hybrid AUC (structured)"], rows,
              "Same folds and training rows; structural features with median imputation and a "
              "missing-structure indicator.")
    t3 = RES / "structure" / "thirdview_cv.csv"
    if t3.exists():
        d = pd.read_csv(t3)
        add_table(doc, "Table S6b. Structure as a third view with the ESM model.",
                  ["Fold", "Sequence ensemble AUC", "Structure ensemble AUC"],
                  [[int(r.fold), f"{r.auc_seq_ensemble:.3f}", f"{r.auc_struct_ensemble:.3f}"]
                   for r in d.itertuples()], "")

    # S7 statistical comparison
    st = json.loads((RES / "stats_vs_published.json").read_text(encoding="utf-8"))
    add_table(doc, "Table S7. Bootstrap comparison with published human performance.",
              ["Quantity", "Value"],
              [["PalmDL human AUC", f"{st['human_auc']:.4f}"],
               ["95% bootstrap CI", f"{st['ci95'][0]:.4f} – {st['ci95'][1]:.4f}"],
               ["n sites / positives", f"{st['n']:,} / {st['n_pos']}"],
               ["P(AUC ≤ 0.900, GPS-Palm human)", f"{st['p_leq_0900']:.2f}"],
               ["P(AUC ≤ 0.855, GPS-Palm general)", f"{st['p_leq_0855']:.3f}"]],
              "4,000 bootstrap replicates over human test sites; the human comparison is a "
              "statistical tie, the general comparison is not.")

    doc.add_heading("Notes on reproducibility", level=2)
    doc.add_paragraph(
        "All supplementary tables are generated from files in results/ by src/make_supplement.py. "
        "The full pipeline (src/01-25) reproduces every number from the UniProt 2026_03 release; "
        "the fold assignments in data/processed/splits_cv.csv are the canonical benchmark folds.")
    doc.save(OUT)
    print("saved", OUT)


if __name__ == "__main__":
    main()
