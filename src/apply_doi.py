"""Fill the Zenodo DOI and repository URL into the manuscript and repository files.

Usage:
    python src/apply_doi.py 10.5281/zenodo.1234567

Edits in place (so the local git repository in dist/github survives):
  manuscript/manuscript.md, dist/github/README.md, dist/github/CITATION.cff,
  dist/github/.zenodo.json
then rebuilds the manuscript DOCX and the Zenodo release zip.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORG = "opus940810"
REPO = "palm-dl"


def main() -> None:
    if len(sys.argv) != 2 or not re.fullmatch(r"10\.5281/zenodo\.\d+", sys.argv[1]):
        raise SystemExit("usage: python src/apply_doi.py 10.5281/zenodo.1234567")
    doi = sys.argv[1]
    edited = []

    ms = ROOT / "manuscript" / "manuscript.md"
    t = ms.read_text(encoding="utf-8")
    if "[ZENODO-DOI]" in t:
        t = t.replace("[ZENODO-DOI]", doi)
        ms.write_text(t, encoding="utf-8")
        edited.append(ms)

    gh = ROOT / "dist" / "github"
    readme = gh / "README.md"
    t = readme.read_text(encoding="utf-8")
    t = t.replace("[ZENODO-DOI]", doi).replace("[ORG]", ORG)
    # keep the archive description accurate: the big files are a release asset,
    # the Zenodo record archives the repository snapshot itself
    t = t.replace(
        "The exact UniProt 2026_03 snapshot, out-of-fold predictions, the proteome-wide\n"
        "score table and the release archive are in the Zenodo record (DOI [ZENODO-DOI])\n"
        "to keep this repository light. Raw data are CC BY 4.0 (UniProt), AlphaFold models\n"
        "CC BY 4.0.",
        f"This repository is archived at Zenodo (DOI {doi}). The exact UniProt 2026_03\n"
        f"snapshot, out-of-fold predictions, the proteome-wide score table and the release\n"
        f"archive are attached to the GitHub release `v1.0.0` as `palmdl-v1.0.0.zip`, because\n"
        f"they are too large for git history. Raw data are CC BY 4.0 (UniProt), AlphaFold\n"
        f"models CC BY 4.0.")
    readme.write_text(t, encoding="utf-8")
    edited.append(readme)

    cit = gh / "CITATION.cff"
    t = cit.read_text(encoding="utf-8")
    t = t.replace("[ZENODO-DOI]", doi).replace("[ORG]", ORG)
    cit.write_text(t, encoding="utf-8")
    edited.append(cit)

    zen = gh / ".zenodo.json"
    data = json.loads(zen.read_text(encoding="utf-8"))
    data.setdefault("related_identifiers", []).append(
        {"relation": "isIdenticalTo", "identifier": f"https://github.com/{ORG}/{REPO}", "scheme": "url"})
    zen.write_text(json.dumps(data, indent=2), encoding="utf-8")
    edited.append(zen)

    print("updated:")
    for p in edited:
        print("  ", p.relative_to(ROOT))
    print(f"\nDOI: {doi}\nrepository: https://github.com/{ORG}/{REPO}")
    print("\nnext:")
    print("  python src/make_docx.py")
    print("  python src/build_release.py")
    print("  cd dist/github && git add -A && git commit -m \"Add Zenodo DOI and repository metadata\"")
    print("  git push")


if __name__ == "__main__":
    main()
