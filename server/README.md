# PalmDL prediction service

Standard-library HTTP service that scores every cysteine of one or more FASTA
proteins with the released PalmDL model (ESM-2 view + motif/physicochemical
random forest, Platt-calibrated). No web framework is required. An optional
structure-aware path uses AlphaFold models when they are available.

## Run

```powershell
$py = "<python.exe>"
$env:PYTHONPATH = "<workspace>\.pylibs"
$env:PALMDL_THREADS = "4"

# interactive server (sequence-only)
& $py palm-dl\server\app.py                 # http://127.0.0.1:8765
& $py palm-dl\server\app.py --port 9000     # custom port

# structure-aware path: AlphaFold models are looked up by accession
& $py palm-dl\server\app.py --structure
& $py palm-dl\server\app.py --structure --no-fetch   # cache only, never download

# batch mode (no server)
& $py palm-dl\server\app.py --batch palm-dl\server\example.fa `
      --out palm-dl\server\example_predictions.csv --threshold 0.5 --structure
```

The first request loads the model (a few seconds); subsequent requests reuse it.

## Endpoints

| Method | Path | Body | Response |
|---|---|---|---|
| GET | `/` | — | HTML form |
| GET | `/health` | — | model metadata, structure availability, cached model count |
| POST | `/predict` | form fields `fasta`, `threshold`, optional `structure`; or JSON `{"fasta": "...", "threshold": 0.5, "structure": true}` | ranked table (HTML) or JSON |

Example:

```powershell
$body = @{ fasta = (Get-Content palm-dl\server\example.fa -Raw); threshold = 0.9; structure = $true } | ConvertTo-Json
Invoke-RestMethod -Uri http://127.0.0.1:8765/predict -Method Post -Body $body -ContentType "application/json"
```

## Structure path

The v1 service exposed an optional structure-aware path backed by a hybrid forest. That model was
trained on the superseded label set, so v2 disables the path: `--structure` is accepted for
backward compatibility but scoring always uses the sequence ensemble. The paired AlphaFold
analysis (forest AUC 0.890 → 0.911) is reported in the manuscript; a retrained hybrid is future
work.

## Choosing a threshold

- `0.36` — natural-prevalence operating point of the deployed logit-average model.
- `0.5` — conservative default for ranking.
- The nested cross-validation benchmark and the deployed model are different algorithms: the
  benchmark selects members per fold and averages within-fold ranks, whereas the service fixes
  RF+ESM and averages logits. Proteome-scan tiers are on a forest-only scale and are not
  comparable with service probabilities.

## Notes

- The nested benchmark chooses members per fold and averages within-fold **ranks**; this service
  fixes RF+ESM and averages **logits**. For the fixed two-view combination the logit average scored
  AUC 0.940 versus 0.944 for the rank average in a sensitivity check; this does not establish that
  the deployed probabilities are calibrated like the benchmark.
- The structure path is retired and the interface no longer exposes it; the AlphaFold analysis is
  an annotated-human paired experiment reported in the manuscript, not a deployed feature.

- Structure helps where annotations already exist (forest AUC 0.890 → 0.910;
  hybrid full model 0.946 → 0.949, manuscript §2.12) but does not rescue the
  motif-atypical literature sites, so the default remains sequence-only.
- Structure lookup needs the record name to be a UniProt accession; arbitrary
  headers will use the sequence path.
- No authentication and no upload limit beyond memory: intended for local or
  intranet use. Put it behind a reverse proxy if exposed.
- Output columns: protein, position, probability, 31-residue window, path, rank.
