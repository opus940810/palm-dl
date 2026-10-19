"""Disable the stale structure path in the v2 service and refresh its README."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

app = ROOT / "server" / "app.py"
t = app.read_text(encoding="utf-8")
old = '''    Handler.use_structure = args.structure
    Handler.allow_fetch = not args.no_fetch'''
new = '''    if args.structure:
        print("NOTE: the structure path is disabled in v2 (the hybrid model was trained on the "
              "superseded labels); scoring uses the sequence ensemble.")
    Handler.use_structure = False
    Handler.allow_fetch = not args.no_fetch'''
assert old in t
app.write_text(t.replace(old, new), encoding="utf-8")

readme = ROOT / "server" / "README.md"
r = readme.read_text(encoding="utf-8")
r = r.replace('''## Choosing a threshold

- `0.34` — natural-prevalence operating point of the sequence model.
- `0.5` — convenient default for ranking; roughly the balanced operating point.
- `0.957` — the proteome-wide **5% FDR** tier (manuscript Table 8); use this for a
  high-confidence shortlist.
- `0.864` — the 30% FDR tier, for discovery rather than confirmation.''',
'''## Choosing a threshold

- `0.36` — natural-prevalence operating point of the released model.
- `0.5` — conservative default for ranking.
- The proteome-scan tiers in the manuscript (Table 8) are on the **forest-view** scale and are
  not comparable with service probabilities; use them only to rank candidates inside that table.''')
r = r.replace('''## Structure-aware path

Enabled with `--structure` (or `"structure": true` in a JSON request). For each
record the service takes the FASTA header's first token as a UniProt accession,
looks for a cached AlphaFold PDB (`results/structure/af/<ACC>.pdb`) and, unless
`--no-fetch` is given, downloads it from the AlphaFold DB API if missing. For
every cysteine whose residue is modelled, the hybrid forest
(`models/rf_hybrid_struct.joblib`, trained on all 71,594 sites with a
missing-structure indicator) replaces the sequence forest, and its logit-average
with the ESM view is calibrated by a dedicated Platt model fitted on out-of-fold
predictions. Cysteines without a modelled residue, and proteins without a model,
fall back to the sequence path; the chosen path is reported per protein in the
`path` column.

Training-time standardisation of the structural block is stored inside the
artifact (`scaler_mean`, `scaler_std`, `median`), so inference uses exactly the
training transformation.''',
'''## Structure path

The v1 service exposed an optional structure-aware path backed by a hybrid forest. That model was
trained on the superseded label set, so v2 disables the path: `--structure` is accepted for
backward compatibility but scoring always uses the sequence ensemble. The paired AlphaFold
analysis (forest AUC 0.890 → 0.911) is reported in the manuscript; a retrained hybrid is future
work.''')
readme.write_text(r, encoding="utf-8")
print("server updated")
