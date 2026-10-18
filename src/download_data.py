"""Download raw data from UniProt REST (streaming endpoint, one request/table).

Outputs (TSV, one row per protein):
  data/raw/uniprot_lipid_reviewed.tsv  - Swiss-Prot entries with >=1 lipidation feature
  data/raw/uniprot_negatives.tsv       - reviewed human entries with no lipidation feature
  data/raw/download_meta.json          - query + UniProt release used
"""
from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

BASE = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = ("accession,sequence,organism_id,protein_name,gene_names,ft_lipid,length,"
          "ft_signal,ft_topo_dom,reviewed")


def stream_tsv(query: str, out_path: Path, timeout: int = 900) -> int:
    params = {"query": query, "format": "tsv", "fields": FIELDS}
    url = f"{BASE}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp, out_path.open("wb") as fh:
        release = resp.headers.get("x-uniprot-release")
        n = 0
        for chunk in iter(lambda: resp.read(1 << 20), b""):
            fh.write(chunk)
            n += chunk.count(b"\n")
    print(f"{out_path.name}: {n - 1} rows, UniProt release {release}", flush=True)
    return n - 1


def main() -> None:
    q_lipid = "(ft_lipid:*) AND (reviewed:true)"
    q_neg = "(reviewed:true) AND (organism_id:9606) AND (NOT ft_lipid:*)"
    n1 = stream_tsv(q_lipid, RAW / "uniprot_lipid_reviewed.tsv")
    n2 = stream_tsv(q_neg, RAW / "uniprot_negatives.tsv")
    meta = {"queries": {"lipid": q_lipid, "negatives": q_neg},
            "rows": {"lipid": n1, "negatives": n2}, "fields": FIELDS,
            "source": "UniProtKB REST stream"}
    (RAW / "download_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("done")


if __name__ == "__main__":
    sys.exit(main())
