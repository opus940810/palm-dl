"""Preserve the git dir when regenerating the bundle, and add a changelog."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

rp = ROOT / "src" / "make_review_package.py"
t = rp.read_text(encoding="utf-8").split("\n")
for i, line in enumerate(t):
    if "71,594" in line:
        t[i] = line.replace("71,594", "67,266")
rp.write_text("\n".join(t), encoding="utf-8")

gh = ROOT / "src" / "make_github_bundle.py"
g = gh.read_text(encoding="utf-8")
old = '''    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)'''
new = '''    if OUT.exists():
        # keep .git so the published history survives regeneration
        for child in OUT.iterdir():
            if child.name == ".git":
                continue
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    OUT.mkdir(parents=True, exist_ok=True)'''
assert old in g
g = g.replace(old, new)

changelog = '''

CHANGELOG = """# Changelog

## v2 (19 October 2026) - label correction

- The v1 dataset filtered lipidation notes with the substring `"palmitoyl"`, admitting
  2,381 N-/O-palmitoyl sites (35% of positives). v2 keeps only S-palmitoyl cysteine on a
  cysteine residue: **4,424 positives** (449 experimental, 3,975 inferred), 67,266 sites,
  21,568 proteins, 786 organisms, 557 human positives.
- Splits: every protein (including background negatives) is clustered; groups sharing an
  exact window are merged. Cross-fold identical windows: **0** (v1: 751 windows).
- Results corrected throughout: cluster CV AUC 0.944 (0.940-0.948), MCC 0.611; human
  0.914 (0.900-0.927); species holdout 0.885; hard negatives 0.956.
- The "UniProt annotation error" claim for the 66 non-cysteine sites is withdrawn: those
  were legitimate N6-/N-/O-palmitoyl annotations admitted by our filter.
- Cross-species matrix reduced to six species; bacterial N-palmitoyl sites removed.
- Proteome scan rescored with the forest view (ESM view not tractable at proteome scale
  on the available CPU host); tiers reported as candidate rankings.
- ESM-2 view omitted from the scan; the released prediction service still uses the full
  sequence ensemble.
"""
'''
g = g.replace("\n\ndef main() -> None:", "\n" + changelog + "\n\ndef main() -> None:")
g = g.replace('(OUT / ".gitignore").write_text(GITIGNORE, encoding="utf-8")',
              '(OUT / "CHANGELOG.md").write_text(CHANGELOG, encoding="utf-8")\n'
              '    (OUT / ".gitignore").write_text(GITIGNORE, encoding="utf-8")')
gh.write_text(g, encoding="utf-8")
print("patched")
