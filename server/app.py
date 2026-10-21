"""PalmDL prediction service (standard library only, no web framework).

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
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import joblib  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

from af_features import fetch as fetch_af  # noqa: E402
from af_features import parse_pdb  # noqa: E402
from config import MODELS, RESULTS, WINDOW  # noqa: E402
from features import context_features, encode_flat  # noqa: E402
from final_model import embed_windows, load, logit  # noqa: E402

AF_DIR = RESULTS / "structure" / "af"
STRUCT_FEATS = ["plddt_site", "plddt_window", "packing"]

SEQ = {}
STRUCT = {}


def get_model():
    if not SEQ:
        SEQ.update(zip(("rf", "net", "cal", "meta"), load()))
    return SEQ


def get_structure_model():
    if not STRUCT:
        obj = joblib.load(MODELS / "rf_hybrid_struct.joblib")
        STRUCT.update(rf=obj["model"], median=obj["median"],
                      mean=obj["scaler_mean"], std=obj["scaler_std"],
                      cal=joblib.load(MODELS / "calibrator_hybrid.joblib"))
    return STRUCT


def parse_fasta(text: str) -> list[tuple[str, str]]:
    records, name, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name = line[1:].split()[0] if len(line) > 1 else f"seq{len(records) + 1}"
            seq = []
        else:
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def windows_for(seq: str) -> list[dict]:
    seq = "".join(c for c in seq.upper() if c.isalpha())
    rows = []
    for i, aa in enumerate(seq):
        if aa != "C":
            continue
        pos = i + 1
        s = pos - 1
        lo, hi = s - 15, s + 16
        left_pad, right_pad = max(0, -lo), max(0, hi - len(seq))
        core = seq[max(0, lo): min(len(seq), hi)]
        rows.append({"pos": pos, "window": "X" * left_pad + core + "X" * right_pad,
                     "protein_length": len(seq), "compartment": "unknown"})
    return rows


def pdb_features(acc: str, positions: list[int]) -> dict[int, dict]:
    """pLDDT and C-alpha packing for requested positions, if modelled."""
    path = AF_DIR / f"{acc}.pdb"
    if not path.exists():
        return {}
    ca = parse_pdb(str(path))
    if not ca:
        return {}
    nums = np.array(sorted(ca))
    xyz = np.array([ca[n][:3] for n in nums])
    pld = np.array([ca[n][3] for n in nums])
    order = {n: i for i, n in enumerate(nums)}
    out = {}
    for pos in positions:
        cand = [n for n in (pos, pos - 1, pos + 1, pos - 2, pos + 2) if n in order]
        if not cand:
            continue
        p = cand[0]
        i = order[p]
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        neigh = (d < 10.0) & (np.abs(np.arange(len(nums)) - i) >= 4)
        win = np.abs(nums - p) <= 5
        out[pos] = {"plddt_site": float(pld[i]), "plddt_window": float(pld[win].mean()),
                    "packing": int(neigh.sum())}
    return out


def combine(p_rf: np.ndarray, p_esm: np.ndarray) -> np.ndarray:
    return 0.5 * (logit(p_rf) + logit(p_esm))


def predict_fasta(text: str, threshold: float, use_structure: bool = False,
                  allow_fetch: bool = True) -> dict:
    records = parse_fasta(text)
    if not records:
        raise ValueError("no FASTA records found")
    seq_model = get_model()
    rf, net, cal = seq_model["rf"], seq_model["net"], seq_model["cal"]
    results = []
    n_struct_sites = 0
    for name, seq in records:
        rows = windows_for(seq)
        if not rows:
            results.append({"protein": name, "length": len(seq), "path": "empty", "sites": []})
            continue
        df = pd.DataFrame(rows)
        X = encode_flat(df.window.tolist())
        C = context_features(df)
        p_rf = rf.predict_proba(np.hstack([X, C]))[:, 1]
        E = embed_windows(df.window.tolist(), bs=192)
        with torch.no_grad():
            p_esm = torch.sigmoid(net(torch.as_tensor(E), torch.as_tensor(C))).numpy()
        z = combine(p_rf, p_esm)
        p = cal.predict_proba(z.reshape(-1, 1))[:, 1]
        path = "sequence"
        if use_structure:
            if allow_fetch and not (AF_DIR / f"{name}.pdb").exists():
                fetch_af(name)
            feats = pdb_features(name, [int(r.pos) for r in df.itertuples()])
            if feats:
                st = get_structure_model()
                S = df[[ "pos"]].copy()
                vals = [feats.get(int(r.pos)) for r in df.itertuples()]
                S = pd.DataFrame([v if v else {f: np.nan for f in STRUCT_FEATS} for v in vals])
                S = S.fillna(pd.Series(st["median"]))
                has = np.array([0.0 if v is None else 1.0 for v in vals]).reshape(-1, 1)
                S = np.column_stack([S.values.astype(float), has])
                S = (S - st["mean"]) / st["std"]
                p_rf_h = st["rf"].predict_proba(np.hstack([X, C, S]))[:, 1]
                z_h = combine(p_rf_h, p_esm)
                p_h = st["cal"].predict_proba(z_h.reshape(-1, 1))[:, 1]
                modelled = np.array([v is not None for v in vals])
                p = np.where(modelled, p_h, p)
                path = "structure" if modelled.any() else "sequence (no AlphaFold model)"
                n_struct_sites += int(modelled.sum())
        hits = []
        for i, r in df.assign(p=p).sort_values("p", ascending=False).iterrows():
            hits.append({"position": int(r.pos), "probability": round(float(r.p), 4),
                         "context": r.window,
                         "passes_threshold": bool(r.p >= threshold)})
        results.append({"protein": name, "length": len(seq), "path": path, "sites": hits})
    flat = [s for r in results for s in r["sites"]]
    flat.sort(key=lambda s: -s["probability"])
    for rank, s in enumerate(flat, 1):
        s["rank"] = rank
    return {"n_proteins": len(records), "n_cysteines": len(flat),
            "threshold": threshold, "structure": bool(use_structure),
            "n_structure_sites": n_struct_sites, "results": results, "model": seq_model["meta"]}


HTML = """<!doctype html><html><head><meta charset="utf-8">
<title>PalmDL — S-palmitoylation site prediction</title>
<style>
body{font-family:system-ui,Arial,sans-serif;max-width:900px;margin:40px auto;padding:0 16px;color:#222}
h1{font-size:22px} textarea{width:100%;height:180px;font-family:ui-monospace,Consolas,monospace;font-size:12px}
input[type=number]{width:90px} .hint{color:#666;font-size:13px} table{border-collapse:collapse;margin-top:16px;width:100%}
th,td{border-bottom:1px solid #ddd;padding:6px 8px;font-size:13px;text-align:left} th{background:#f6f6f6}
code{background:#f3f3f3;padding:1px 4px;border-radius:3px} .hi{color:#b22222;font-weight:600}
</style></head><body>
<h1>PalmDL — protein S-palmitoylation site prediction</h1>
<p class="hint">Paste one or more protein sequences in FASTA format. Every cysteine is scored and
the table is ranked by probability. The service uses the fixed RF+ESM logit-average model with a
Platt calibration; it is not the nested cross-validation benchmark. The natural-prevalence
operating point is 0.36; higher thresholds are more conservative. Proteome-scan tiers in the
manuscript are on a different (forest-only) scale and are not comparable with these probabilities.</p>
<form method="post" action="/predict">
<textarea name="fasta" placeholder=">P25101 endothelin receptor A\nMN..."></textarea><br><br>
<label>threshold <input type="number" name="threshold" step="0.01" value="0.5"></label>
&nbsp;<button type="submit">Predict</button>
</form></body></html>"""


def render(result: dict) -> str:
    rows = []
    for rec in result["results"]:
        for s in rec["sites"]:
            if not s["passes_threshold"]:
                continue
            rows.append(f"<tr><td>{rec['protein']}</td><td>C{s['position']}</td>"
                        f"<td>{s['probability']:.3f}</td><td><code>{s['context']}</code></td>"
                        f"<td>{rec['path']}</td></tr>")
    body = "".join(rows) or "<tr><td colspan=5>No cysteine passed the threshold.</td></tr>"
    return HTML.replace("</form></body></html>", f"""</form>
<p class="hint">{result['n_proteins']} proteins, {result['n_cysteines']} cysteines;
threshold {result['threshold']}; structure={result['structure']}
({result['n_structure_sites']} cysteines with structural features); {len(rows)} passing.</p>
<table><tr><th>protein</th><th>site</th><th>P</th><th>window</th><th>path</th></tr>
{body}</table></body></html>""")


class Handler(BaseHTTPRequestHandler):
    use_structure = False
    allow_fetch = True

    def _send(self, code: int, body: str, ctype: str = "text/html; charset=utf-8") -> None:
        data = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:  # noqa: N802
        if self.path.startswith("/health"):
            meta = get_model()["meta"]
            payload = {"status": "ok", "model": meta,
                       "structure_path_available": self.use_structure,
                       "cached_alphafold_models": len(list(AF_DIR.glob("*.pdb")))}
            self._send(200, json.dumps(payload), "application/json")
        elif self.path in ("/", "/index.html"):
            self._send(200, HTML)
        else:
            self._send(404, "not found")

    def do_POST(self) -> None:  # noqa: N802
        if not self.path.startswith("/predict"):
            self._send(404, "not found")
            return
        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length).decode("utf-8", "ignore")
        ctype = self.headers.get("Content-Type", "")
        try:
            if "application/json" in ctype:
                payload = json.loads(raw)
                fasta = payload.get("fasta", "")
                threshold = float(payload.get("threshold", 0.5))
                use_struct = bool(payload.get("structure", self.use_structure))
            else:
                form = parse_qs(raw)
                fasta = form.get("fasta", [""])[0]
                threshold = float(form.get("threshold", ["0.5"])[0])
                use_struct = bool(form.get("structure", [""])[0]) or self.use_structure
            # the legacy structure path is retired: the hybrid artifact was trained on
            # superseded labels, so requests can never re-enable it
            result = predict_fasta(fasta, threshold, False, self.allow_fetch)
        except Exception as exc:  # noqa: BLE001
            self._send(400, f"error: {exc}")
            return
        if "application/json" in ctype:
            self._send(200, json.dumps(result), "application/json")
        else:
            self._send(200, render(result))


def batch(fasta_path: Path, out_path: Path, threshold: float, use_structure: bool,
          allow_fetch: bool) -> None:
    result = predict_fasta(fasta_path.read_text(encoding="utf-8"), threshold,
                           use_structure, allow_fetch)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["protein", "position", "probability", "window", "path", "rank"])
        for rec in result["results"]:
            for s in rec["sites"]:
                w.writerow([rec["protein"], s["position"], s["probability"], s["context"],
                            rec["path"], s["rank"]])
    print(f"{sum(len(r['sites']) for r in result['results'])} cysteines -> {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--structure", action="store_true", help="enable the structure-aware path")
    ap.add_argument("--no-fetch", action="store_true", help="never download AlphaFold models")
    ap.add_argument("--batch")
    ap.add_argument("--out", default="predictions.csv")
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()
    if args.structure:
        print("NOTE: the structure path is disabled in v2 (the hybrid model was trained on "
              "the superseded labels); scoring uses the sequence ensemble.")
        args.structure = False
    Handler.use_structure = args.structure
    Handler.allow_fetch = not args.no_fetch
    if args.batch:
        batch(Path(args.batch), Path(args.out), args.threshold, args.structure, not args.no_fetch)
        return
    print(f"PalmDL server on http://{args.host}:{args.port} "
          f"(structure={'on' if args.structure else 'off'}, fetch={'off' if args.no_fetch else 'on'})")
    ThreadingHTTPServer((args.host, args.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
