"""Build the supplementary document (Tables S1-S7) for the S-only rebuild."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from docx import Document
from docx.shared import Pt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "manuscript" / "PalmDL_supplementary.docx"


def add_table(doc, title: str, header: list[str], rows: list[list], note: str = "", size: float = 7.0):
    doc.add_heading(title, level=2)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        r = c.paragraphs[0].add_run(h)
        r.bold = True
        r.font.size = Pt(size)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row[: len(header)]):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(str(v))
            run.font.size = Pt(size)
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
    doc.add_paragraph("This supplement documents the label correction, per-fold metrics, the leakage "
                      "audit, the curated literature set, the proteome scan tiers, the structure "
                      "comparison and the statistical comparison with published values.")

    add_table(doc, "Table S1. Composition of the first label set and the correction.",
              ["Annotation chemistry", "Unique sites", "Share", "Kept in v2"],
              [["S-palmitoyl cysteine", "4,424", "65.0%", "yes"],
               ["N-palmitoyl cysteine", "2,315", "34.0%", "no"],
               ["N6-palmitoyl lysine", "34", "0.50%", "no"],
               ["N-palmitoyl glycine", "18", "0.26%", "no"],
               ["O-palmitoyl serine", "9", "0.13%", "no"],
               ["O-palmitoyl threonine", "5", "0.07%", "no"],
               ["Total admitted by the first filter", "6,805", "100%", "4,424 kept"]],
              "The first filter matched any lipidation note containing 'palmitoyl'. The corrected "
              "filter requires an S-palmitoyl cysteine note on a cysteine residue.")

    frames = []
    for f in ("metrics_cv_base.csv", "metrics_cv_cnn.csv", "metrics_cv_esm.csv",
              "metrics_cv_ensemble.csv"):
        p = RES / f
        if p.exists():
            d = pd.read_csv(p)
            frames.append(d[d.protocol == "test_natural"])
    if frames:
        m = pd.concat(frames)
        agg = (m.groupby("model")[["mcc", "auc", "aupr", "acc", "sn", "sp"]]
               .agg(["mean", "std"]).round(4))
        rows = [[model] + [f"{r[(c, 'mean')]:.3f} ± {r[(c, 'std')]:.3f}"
                           for c in ["auc", "aupr", "mcc", "acc", "sn", "sp"]]
                for model, r in agg.iterrows()]
        add_table(doc, "Table S2. Per-fold cross-validation metrics (mean ± SD over five folds).",
                  ["Model", "AUC", "AUPR", "MCC", "Accuracy", "Sn", "Sp"], rows,
                  "Cluster-held-out folds, all species, natural prevalence; thresholds chosen on "
                  "inner validation folds.")

    s3 = pd.read_csv(RES / "leakage_ladder.csv")
    add_table(doc, "Table S3. Split-protocol audit (same random forest and features).",
              ["Protocol", "Train pos.", "Test sites", "Test pos.", "AUC", "MCC", "Balanced MCC",
               "Shared windows"],
              [[r.protocol, int(r.n_train_pos), int(r.n_test), int(r.n_test_pos), f"{r.auc:.3f}",
                f"{r.mcc_nat:.3f}", f"{r.mcc_bal:.3f}", f"{r.leaked_pos_windows * 100:.1f}%"]
               for r in s3.itertuples()],
              "'Shared windows' is the fraction of test positives whose exact 31-residue window "
              "occurs among training positives; cluster and species splits share none.")

    cur = pd.read_csv(RES / "literature" / "curated_sites.csv")
    cur = cur[cur.curated]
    names = protein_names()
    rows = [[i, r.acc, f"C{r.pos}", r.token, str(names.get(r.acc, ""))[:60], r.pmcid]
            for i, r in enumerate(cur.itertuples(), 1)]
    add_table(doc, "Table S4. Curated literature sites (99 of 136 extracted mentions).",
              ["#", "Accession", "Site", "Token", "Protein", "PMCID"], rows,
              "Accepted when the sentence attributes palmitoylation to that protein at that "
              "cysteine, or a cysteine mutant loses palmitoylation; rejected otherwise.",
              size=6.0)

    fdr = pd.read_csv(RES / "scan" / "fdr_summary.csv")
    add_table(doc, "Table S5a. Proteome-wide scan tiers (forest view of the released model).",
              ["Tier", "Sites", "Proteins", "Known recovered", "Recall", "Novel", "Threshold"],
              [[f"{r.fdr:.0%}", int(r.n_predictions), int(r.n_proteins), int(r.known_recovered),
                f"{r.known_recall:.1%}", int(r.novel), f"{r.threshold:.3f}"] for r in fdr.itertuples()],
              "262,286 cysteines in 19,814 reviewed human proteins; 20,000 shuffled-flank decoys. "
              "The decoy tail is thin, so tiers are candidate rankings rather than proven error rates.")
    ce = pd.read_csv(RES / "scan" / "compartment_enrichment.csv")
    add_table(doc, "Table S5b. Topology enrichment of the 5% tier.",
              ["Compartment", "Share in shortlist", "Share in all cysteines", "Enrichment"],
              [[r.compartment, f"{r.share_shortlist:.1%}", f"{r.share_all:.1%}", f"{r.enrichment:.1f}×"]
               for r in ce.itertuples()],
              "Topology from UniProt; 'unknown' dominates the background and is depleted.")

    pc = pd.read_csv(RES / "structure" / "paired_cv.csv")
    rows = []
    for fold, g in pc.groupby("fold"):
        b = g[g.model == "base"].iloc[0]
        s = g[g.model == "base+structure"].iloc[0]
        rows.append([int(fold), int(b.n), int(b.n_pos), f"{b.auc:.3f}", f"{s.auc:.3f}",
                     f"{s.auc - b.auc:+.3f}"])
    add_table(doc, "Table S6. AlphaFold features on human sites, paired per fold.",
              ["Fold", "Sites", "Positives", "Base AUC", "+ structure AUC", "Δ"], rows,
              "Same random forest and cluster-held-out folds with and without pLDDT and C-alpha "
              "packing features; mean Δ +0.021, positive in four of five folds.")

    st = json.loads((RES / "stats_vs_published.json").read_text(encoding="utf-8"))
    add_table(doc, "Table S7. Bootstrap comparison with published human performance.",
              ["Quantity", "Value"],
              [["PalmDL human AUC (fold mean)", f"{st['human_auc_fold_mean']:.4f}"],
               ["95% within-fold bootstrap CI", f"{st['ci95'][0]:.4f} – {st['ci95'][1]:.4f}"],
               ["n sites / positives", f"{st['n']:,} / {st['n_pos']}"],
               ["P(AUC ≤ 0.900, GPS-Palm human)", f"{st['p_leq_0900']:.3f}"],
               ["P(AUC ≤ 0.855, GPS-Palm general)", f"{st['p_leq_0855']:.3f}"]],
              "4,000 within-fold bootstrap replicates over human test sites; the human comparison "
              "is marginal and its interval touches 0.900.")

    doc.add_heading("Notes on reproducibility", level=2)
    doc.add_paragraph(
        "All supplementary tables are generated from files in results/ by src/make_supplement.py. "
        "The full pipeline reproduces every number from the UniProt 2026_03 release; the fold "
        "assignments in data/processed/splits_cv.csv are the canonical benchmark folds.")
    doc.save(OUT)
    print("saved", OUT)


if __name__ == "__main__":
    main()
