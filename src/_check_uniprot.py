import json
import sys
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")
u = "https://api.crossref.org/works?" + urllib.parse.urlencode(
    {"query.title": "UniProt universal protein knowledgebase", "rows": 6,
     "filter": "from-pub-date:2024-06-01", "select": "title,container-title,issued,volume,page,DOI"})
d = json.load(urllib.request.urlopen(u, timeout=40))
for it in d["message"]["items"]:
    print((it.get("title") or ["?"])[0][:90], "|",
          (it.get("container-title") or ["?"])[0][:35], "|",
          it.get("issued", {}).get("date-parts", [["?"]])[0][0], "|",
          it.get("volume"), it.get("page"), "|", it["DOI"])
