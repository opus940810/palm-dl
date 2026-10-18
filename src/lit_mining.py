"""Literature-mined external test set for PalmDL.

Stage 1  search Europe PMC for recent open-access palmitoylation papers.
Stage 2  fetch full-text XML for a bounded sample.
Stage 3  extract candidate (protein token, cysteine position) mentions from
         sentences containing "palmitoyl", resolving each token to a reviewed
         human UniProt entry and requiring the residue at that position to be C.
Stage 4  score the verified sites with the deployable PalmDL model and compare
         with same-protein cysteines (within-protein ranking).

Only sites that are absent from the UniProt 2026_03 training table count as
"novel"; the rest measure how well the model reproduces curated knowledge.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "literature"
OUT.mkdir(parents=True, exist_ok=True)
H = {"User-Agent": "Mozilla/5.0 (academic literature mining)"}
BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"

STOP = {"the", "this", "these", "those", "with", "from", "that", "were", "was", "and", "for",
        "cys", "cysteine", "palmitoylation", "palmitoylated", "protein", "proteins", "site",
        "sites", "residue", "residues", "mice", "mouse", "human", "cells", "cell", "wild",
        "type", "mutant", "mutants", "figure", "table", "supplementary", "data", "using",
        "shown", "observed", "reported", "previously", "respectively", "both", "also"}

SITE_RE = re.compile(r"(?<![A-Za-z0-9\-])(?:Cys(?:teine)?[\s-]?|C)(\d{1,4})(?!\d)", re.I)
# gene-like token immediately before "Cys<number>" / "C<number>"
LEFT_GENE = re.compile(r"([A-Za-z][A-Za-z0-9\-]{1,12})\s*(?:protein|residue)?\s*[,:]?\s*(?:at\s+)?$")
RIGHT_GENE = re.compile(r"^\s*(?:of|in|on|from)\s+([A-Za-z][A-Za-z0-9\-]{1,12})")


def _get_json(url: str, retries: int = 4):
    for i in range(retries):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=120))
        except Exception as exc:  # noqa: BLE001
            print(f"  retry {i} ({repr(exc)[:60]})", flush=True)
            time.sleep(3 * (i + 1))
    raise RuntimeError("Europe PMC unavailable")


def search(limit: int) -> list[dict]:
    q = ('(palmitoylation OR "S-palmitoylation" OR palmitoylated) AND OPEN_ACCESS:Y '
         "AND FIRST_PDATE:[2020-01-01 TO 2026-12-31]")
    out, cursor = [], "*"
    while len(out) < limit:
        params = {"query": q, "format": "json", "pageSize": 1000, "resultType": "lite",
                  "cursorMark": cursor}
        url = f"{BASE}/search?{urllib.parse.urlencode(params)}"
        d = _get_json(url)
        batch = [r for r in d["resultList"]["result"] if r.get("pmcid")]
        out.extend(batch)
        nxt = d.get("nextCursorMark")
        print(f"  search {len(out)} (hits {d['hitCount']})", flush=True)
        if not nxt or nxt == cursor or not batch:
            break
        cursor = nxt
    seen, uniq = set(), []
    for r in out:
        if r["pmcid"] not in seen:
            seen.add(r["pmcid"])
            uniq.append(r)
    return uniq[:limit]


def fetch_one(pmcid: str) -> str | None:
    path = OUT / "xml" / f"{pmcid}.xml"
    if path.exists():
        return path.read_text(encoding="utf-8", errors="ignore")
    url = f"{BASE}/{pmcid}/fullTextXML"
    try:
        txt = urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=60).read().decode(
            "utf-8", "ignore")
        path.parent.mkdir(exist_ok=True)
        path.write_text(txt, encoding="utf-8")
        return txt
    except Exception:
        return None


def sentences(text: str) -> list[str]:
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"\s+", " ", text)
    return re.split(r"(?<=[.!?])\s+", text)


def extract_texts(texts: list[str], pmcids: list[str]) -> pd.DataFrame:
    rows = []
    for pmcid, text in zip(pmcids, texts):
        if not text:
            continue
        for sent in sentences(text):
            low = sent.lower()
            if "palmitoyl" not in low:
                continue
            for m in SITE_RE.finditer(sent):
                pos = int(m.group(1))
                if not (1 <= pos <= 3000):
                    continue
                left = sent[max(0, m.start() - 40):m.start()]
                right = sent[m.end():m.end() + 40]
                cand = []
                lg = LEFT_GENE.search(left)
                if lg:
                    cand.append(lg.group(1))
                rg = RIGHT_GENE.search(right)
                if rg:
                    cand.append(rg.group(1))
                for token in cand:
                    rows.append({"pmcid": pmcid, "pos": pos, "token": token,
                                 "match": m.group(0), "sentence": sent[:300]})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["token_l"] = df.token.str.lower()
    df = df[~df.token_l.isin(STOP)]
    df = df[df.token.str.len() >= 3]

    def gene_like(t: str) -> bool:
        if re.match(r"^(?:Cys(?:teine)?|C)[- ]?\d+$", t, re.I):
            return False  # a site reference, not a gene
        return bool(re.search(r"\d", t)) or t.isupper() or bool(
            re.match(r"^[a-z]{1,4}[A-Z][A-Za-z0-9]*$", t))

    df = df[df.token.map(gene_like)]
    return df.drop_duplicates(["pmcid", "pos", "token"])


def map_genes(tokens: list[str]) -> dict[str, str]:
    """Map gene/protein tokens to a unique reviewed accession (human preferred)."""
    cache_path = OUT / "gene_map.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    todo = [t for t in dict.fromkeys(tokens) if t not in cache]
    print(f"mapping {len(todo)} new gene tokens ({len(cache)} cached)", flush=True)
    for i in range(0, len(todo), 40):
        chunk = todo[i:i + 40]
        query = " OR ".join(f"(gene:{t} OR protein_name:{t})" for t in chunk)
        params = {"query": f"({query}) AND reviewed:true AND "
                           "(organism_id:9606 OR organism_id:10090 OR organism_id:10116)",
                  "format": "tsv", "size": 500,
                  "fields": "accession,gene_names,organism_id,protein_name"}
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
        txt = None
        for attempt in range(3):
            try:
                txt = urllib.request.urlopen(urllib.request.Request(url, headers=H),
                                             timeout=180).read().decode()
                break
            except Exception as exc:  # noqa: BLE001
                print("  map retry", chunk[:2], repr(exc)[:60])
                time.sleep(2 * (attempt + 1))
        hits: dict[str, list[tuple[str, str]]] = {}
        if txt:
            for line in txt.splitlines()[1:]:
                f = line.split("\t")
                if len(f) < 4:
                    continue
                acc, genes, org, pname = f[0], f[1], f[2], f[3]
                for t in chunk:
                    pat = re.compile(rf"\b{re.escape(t)}\b", re.I)
                    if pat.search(genes) or pat.search(pname):
                        hits.setdefault(t, []).append((acc, org))
        for t in chunk:
            hl = hits.get(t, [])
            human = [a for a, o in hl if o.startswith("9606")]
            if len(human) == 1:
                cache[t] = human[0]
            elif len(hl) == 1:
                cache[t] = hl[0][0]
            else:
                cache[t] = ""
        time.sleep(0.2)
    cache_path.write_text(json.dumps(cache), encoding="utf-8")
    return cache


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=1500)
    ap.add_argument("--stage", default="all")
    args = ap.parse_args()

    meta_path = OUT / "search.json"
    if args.stage in ("all", "search"):
        res = search(args.limit)
        meta_path.write_text(json.dumps(res, indent=1), encoding="utf-8")
        print("saved", len(res), "records")

    recs = json.loads(meta_path.read_text())
    pmcids = [r["pmcid"] for r in recs]
    if args.stage in ("all", "fetch"):
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=8) as ex:
            texts = list(ex.map(fetch_one, pmcids))
        ok = sum(1 for t in texts if t)
        print(f"fetched {ok}/{len(texts)} full texts in {time.time() - t0:.0f}s", flush=True)
    else:
        texts = [fetch_one(p) for p in pmcids]

    cand = extract_texts(texts, pmcids)
    cand.to_csv(OUT / "candidates.csv", index=False)
    print("candidate mentions:", len(cand), "unique tokens:", cand.token.nunique())
    print(cand.head(8).to_string())


if __name__ == "__main__":
    main()
