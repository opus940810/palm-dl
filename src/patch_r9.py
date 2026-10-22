"""Cite the concept DOI only (versions are pinned by the GitHub tag), bump archive."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")
t = t.replace(
    "archived at Zenodo under the concept DOI 10.5281/zenodo.23455872 (version 2.1.0: 10.5281/zenodo.23537478) and available at GitHub tag v2.1",
    "archived at Zenodo under the concept DOI 10.5281/zenodo.23455872 (latest version 2.1.1; the exact snapshot is pinned by GitHub tag v2.1) and available at GitHub tag v2.1")
md.write_text(t, encoding="utf-8")

wp = ROOT / "src" / "write_reproduce.py"
w = wp.read_text(encoding="utf-8")
w = w.replace("Archived at Zenodo under the concept DOI **{CONCEPT}** (version 2.1.0: {VERSION}).",
              "Archived at Zenodo under the concept DOI **{CONCEPT}**; the exact snapshot is pinned by tag v2.1.")
wp.write_text(w, encoding="utf-8")

br = ROOT / "src" / "build_release.py"
b = br.read_text(encoding="utf-8").replace('VERSION = "v2.1.1"', 'VERSION = "v2.1.2"')
br.write_text(b, encoding="utf-8")
print("concept-only citation applied; archive v2.1.2")
