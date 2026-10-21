"""Point every DOI reference at the v2.1 concept DOI and refresh the archive wording."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONCEPT = "10.5281/zenodo.23455872"
VERSION = "10.5281/zenodo.23537476"

# manuscript availability paragraph
md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")
old_start = "The corrected v2.1 code, processed data"
i = t.find(old_start)
if i != -1:
    j = t.find("\n", i)
    new = (f"The corrected v2.1 code, processed data (including the canonical fold assignments), "
           f"out-of-fold predictions and proteome-wide scores are archived at Zenodo under the concept "
           f"DOI {CONCEPT} (version 2.1.0: {VERSION}) and available at GitHub tag v2.1 "
           f"(https://github.com/opus940810/palm-dl). The superseded v1.0.0 release "
           f"(DOI 10.5281/zenodo.23455873) predates the label correction and must not be used for the "
           f"results reported here. The prediction service is available at [SERVICE-URL]. Supplementary "
           f"Tables S1–S7 accompany the manuscript.")
    t = t[:i] + new + t[j:]
md.write_text(t, encoding="utf-8")

# root README / Chinese notes: replace the old version DOI with the concept DOI
for name in ("README.md", "中文说明.md", "docs/REPRODUCE.md"):
    p = ROOT / name
    if not p.exists():
        continue
    s = p.read_text(encoding="utf-8")
    s = s.replace("10.5281/zenodo.23455873", CONCEPT)
    if name == "docs/REPRODUCE.md":
        s = s.replace("The repository tag is **v2.1**;",
                      f"Archived at Zenodo under the concept DOI **{CONCEPT}** (version 2.1.0: {VERSION}).\n"
                      "The repository tag is **v2.1**;")
    p.write_text(s, encoding="utf-8")

# archive version bump
br = ROOT / "src" / "build_release.py"
b = br.read_text(encoding="utf-8").replace('VERSION = "v2.1.0"', 'VERSION = "v2.1.1"')
br.write_text(b, encoding="utf-8")

# bundle template: use the concept DOI in the data/archive paragraph
gb = ROOT / "src" / "make_github_bundle.py"
g = gb.read_text(encoding="utf-8")
g = g.replace("The corrected v2.1 code and processed data live in this repository (tag v2.1). The\n"
              "earlier v1.0.0 Zenodo record predates the label correction and must not be used for the\n"
              "current results; a corrected Zenodo version accompanies the manuscript.",
              f"The corrected v2.1 code and processed data live in this repository (tag v2.1) and are\n"
              f"archived at Zenodo under the concept DOI {CONCEPT} (version 2.1.0: {VERSION}). The earlier\n"
              f"v1.0.0 record (10.5281/zenodo.23455873) predates the label correction and must not be used\n"
              f"for the current results.")
gb.write_text(g, encoding="utf-8")
print("DOI references updated to", CONCEPT)
