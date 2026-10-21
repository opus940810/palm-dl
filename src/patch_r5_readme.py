"""Fix the remaining structure-promise block in the root README."""
from pathlib import Path

p = Path(__file__).resolve().parents[1] / "README.md"
t = p.read_text(encoding="utf-8")
old = "`--structure` enables the structure-aware path: records whose header is a UniProt accession are\nlooked up in the AlphaFold cache (or downloaded, unless `--no-fetch`), modelled cysteines are scored\nwith the hybrid forest and a dedicated out-of-fold calibrator, and everything else falls back to the\nsequence path; the chosen path is reported per protein. Verified end to end (P25101 -> structure,\nsynthetic record -> sequence fallback)."
new = ("The service is sequence-only. The legacy structure-aware path is retired in both the CLI and HTTP\n"
       "handlers (the hybrid forest was trained on the superseded labels); `--structure` is accepted for\n"
       "backward compatibility but ignored. The AlphaFold analysis in the paper is a paired experiment on\n"
       "annotated human sites, not a deployed feature.")
if old not in t:
    for i, line in enumerate(t.split("\n")):
        if line.startswith("`--structure` enables the structure-aware path"):
            print("found at line", i)
    raise SystemExit("block not found")
t = t.replace(old, new)
t = t.replace("## Structure features and the structure-aware prototype", "## Structure analysis on annotated human sites")
p.write_text(t, encoding="utf-8")
print("README structure block replaced")
