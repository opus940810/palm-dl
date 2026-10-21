"""Remove the last FDR-controlled wording and clean the service docstring."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

md = ROOT / "manuscript" / "manuscript.md"
t = md.read_text(encoding="utf-8")
if "FDR-controlled" in t:
    for line in t.split("\n"):
        if "FDR-controlled" in line:
            print("context:", line[:160])
    t = t.replace("FDR-controlled", "validation-target")
    md.write_text(t, encoding="utf-8")
    print("replaced")

app = ROOT / "server" / "app.py"
s = app.read_text(encoding="utf-8")
start = s.index('"""PalmDL prediction service')
end = s.index('"""', start + 3) + 3
new = '''"""PalmDL prediction service (standard library only, no web framework).

Scoring uses the fixed RF+ESM logit-average model with a Platt calibration. This is
not the nested cross-validation benchmark: the benchmark chooses ensemble members
per fold and averages within-fold ranks, whereas the service fixes RF+ESM and
averages logits. The legacy structure-aware path is retired; requests cannot
enable it.

Run:
    python app.py                          # http://127.0.0.1:8765
    python app.py --batch in.fa --out sites.csv

Endpoints:
    GET  /            HTML form
    POST /predict     form or JSON {"fasta": "...", "threshold": 0.5}
    GET  /health      model metadata
"""'''
app.write_text(s[:start] + new + s[end:], encoding="utf-8")
print("docstring updated")
