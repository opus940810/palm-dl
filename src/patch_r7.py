"""Seventh-round fixes: version DOI, raw snapshot, scan view, metadata templates."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONCEPT = "10.5281/zenodo.23455872"
VERSION = "10.5281/zenodo.23537478"

# 1) manuscript version DOI
md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8").replace("10.5281/zenodo.23537476", VERSION)
md.write_text(t, encoding="utf-8")

# 2) scan default view
ss = ROOT / "src" / "scan_score.py"
s = ss.read_text(encoding="utf-8")
s = s.replace('ap.add_argument("--view", default="esm", choices=["esm", "rf"],',
              'ap.add_argument("--view", default="rf", choices=["esm", "rf"],')
ss.write_text(s, encoding="utf-8")

# 3) release includes the raw snapshot and the manual curation sheet
br = ROOT / "src" / "build_release.py"
b = br.read_text(encoding="utf-8")
if '"data/raw/uniprot_lipid_reviewed.tsv"' not in b:
    b = b.replace('"data/processed", "data/raw/download_meta.json",',
                  '"data/processed", "data/raw", "results/literature/review_sheet.tsv",')
br.write_text(b, encoding="utf-8")

# 4) reproduce doc: version DOI, ESM checkpoint, scan view
wp = ROOT / "src" / "write_reproduce.py"
w = wp.read_text(encoding="utf-8")
w = w.replace('python src/scan_score.py\n', 'python src/scan_score.py --view rf   # forest view used in the paper\n')
w = w.replace('The repository tag is **v2.1**;',
              f"Archived at Zenodo under the concept DOI **{CONCEPT}** (version 2.1.0: {VERSION}).\n"
              "The repository tag is **v2.1**;")
w = w.replace("## Release identity (v2.1)",
              """## ESM-2 checkpoint

The language-model view needs the official fair-esm weights. Download once with
`python -c "import esm; esm.pretrained.esm2_t12_35M_UR50D()"` and place the checkpoint at
`palm-dl/.cache/esm2_t12_35M_UR50D.pt` (the path `src/esm_embed.py` expects). No other
model weights are required; `src/af_features.py` downloads AlphaFold models on demand.

## Release identity (v2.1)""")
wp.write_text(w, encoding="utf-8")

# 5) bundle templates: real author, concept DOI relation, snapshot note
gb = ROOT / "src" / "make_github_bundle.py"
g = gb.read_text(encoding="utf-8")
g = g.replace('"creators": [{"name": "[Family name], [Given name]", "affiliation": "[Affiliation]"}],',
              '"creators": [{"name": "Ruan, Yun", "affiliation": "[Affiliation]"}],')
g = g.replace('"access_right": "open",',
              '"access_right": "open",\n    "related_identifiers": [{"relation": "isSupplementTo", "identifier": "https://github.com/opus940810/palm-dl/tree/v2.1", "scheme": "url"}],')
g = g.replace('  - family-names: "[to be completed]"\n    given-names: "[to be completed]"',
              '  - family-names: "Ruan"\n    given-names: "Yun"')
g = g.replace('The corrected v2.1 code and processed data live in this repository (tag v2.1) and are',
              'The exact UniProt 2026_03 snapshot used for the paper is included under `data/raw/`; re-running\n'
              '`src/download_data.py` refreshes it from the live API, whose release may differ. The corrected\n'
              'v2.1 code and processed data live in this repository (tag v2.1) and are')
gb.write_text(g, encoding="utf-8")
print("seventh-round patch applied")
