# PalmDL prediction service

Standard-library HTTP service that scores every cysteine of one or more FASTA proteins with the
released PalmDL model (ESM-2 view + motif/physicochemical random forest, Platt-calibrated). The
service is **sequence-only**; the legacy structure path is retired.

## Run

```powershell
$py = "<python.exe>"
$env:PYTHONPATH = "<workspace>\.pylibs"
$env:PALMDL_THREADS = "4"

# interactive server
& $py palm-dl\server\app.py                 # http://127.0.0.1:8765
& $py palm-dl\server\app.py --port 9000     # custom port

# batch mode (no server)
& $py palm-dl\server\app.py --batch palm-dl\server\example.fa `
      --out palm-dl\server\example_predictions.csv --threshold 0.36
```

The first request loads the model (a few seconds); subsequent requests reuse it.

## Endpoints

| Method | Path | Body | Response |
|---|---|---|---|
| GET | `/` | — | HTML form |
| GET | `/health` | — | model metadata |
| POST | `/predict` | form fields `fasta`, `threshold`; or JSON `{"fasta": "...", "threshold": 0.36}` | ranked table (HTML) or JSON |

Example:

```powershell
$body = @{ fasta = (Get-Content palm-dl\server\example.fa -Raw); threshold = 0.36 } | ConvertTo-Json
Invoke-RestMethod -Uri http://127.0.0.1:8765/predict -Method Post -Body $body -ContentType "application/json"
```

## Structure

The v1 service exposed an optional structure-aware path backed by a hybrid forest trained on the
superseded label set. v2 closes that path in both the CLI and HTTP handlers; `--structure` is
accepted for backward compatibility but ignored. The AlphaFold analysis in the paper is a paired
experiment on annotated human sites (forest AUC 0.890 → 0.911), not a deployed feature.

## Choosing a threshold

- `0.36` — natural-prevalence operating point of the deployed logit-average model.
- `0.5` — conservative default for ranking.
- The nested cross-validation benchmark and the deployed model are different algorithms: the
  benchmark selects members per fold and averages within-fold ranks, whereas the service fixes
  RF+ESM and averages logits. Proteome-scan tiers are on a forest-only scale and are not
  comparable with service probabilities.

## Notes

- For the fixed two-view combination the logit average scored AUC 0.940 versus 0.944 for the rank
  average in a sensitivity check (difference about 0.004); this does not establish that deployed
  probabilities are calibrated like the benchmark, whose primary estimate is the nested 0.938.
- No authentication and no upload limit beyond memory: intended for local or intranet use. Put it
  behind a reverse proxy if exposed.
