import json
import sys
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")


def q(s):
    u = "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.bibliographic": s, "rows": 2,
         "select": "title,author,container-title,issued,volume,page,DOI"})
    d = json.load(urllib.request.urlopen(u, timeout=40))
    for it in d["message"]["items"]:
        au = it.get("author", [{}])[0].get("family", "?")
        t = (it.get("title") or ["?"])[0][:95]
        c = (it.get("container-title") or ["?"])[0][:45]
        print(" |", au, "|", t, "|", c, "|",
              it.get("issued", {}).get("date-parts", [["?"]])[0][0], "|",
              it.get("volume"), it.get("page"), "|", it["DOI"])


for s in ["UniProt the Universal Protein Knowledgebase in 2025",
          "NBA-Palm prediction of palmitoylation site implemented in Naive Bayes algorithm",
          "SwissPalm protein palmitoylation database",
          "CysModDB comprehensive database of cysteine modifications"]:
    print("==", s)
    q(s)
