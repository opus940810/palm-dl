# Hosting the PalmDL service

The service is a single standard-library Python process; it needs the model
artifacts (`models/`) and, for the structure path, an AlphaFold cache
(`results/structure/af/`). Both are excluded from the Zenodo archive for size and
are regenerated in minutes.

## Option A — a VM with systemd (simplest)

```bash
# on the server
git clone <repository> palm-dl && cd palm-dl
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python src/download_data.py && python src/build_dataset.py && python src/splits.py
python src/esm_embed.py          # ~20 min CPU, needs the ESM-2 35M weights
python src/final_model.py        # ~3 min
python src/structure_hybrid.py   # optional, ~10 min, for the structure path
```

`/etc/systemd/system/palmdl.service`:

```ini
[Unit]
Description=PalmDL prediction service
After=network.target

[Service]
WorkingDirectory=/opt/palm-dl
Environment=PALMDL_THREADS=4
ExecStart=/opt/palm-dl/.venv/bin/python server/app.py --host 127.0.0.1 --port 8765 --structure
Restart=always
User=palmdl

[Install]
WantedBy=multi-user.target
```

Put `deploy/nginx.conf` in `/etc/nginx/sites-available/palmdl`, adjust
`server_name`, then `systemctl enable --now palmdl && nginx -s reload`.
HTTPS: `certbot --nginx -d your.domain`.

## Option B — Docker

```bash
cd palm-dl/deploy
docker compose up -d --build
curl http://localhost:8765/health
```

`docker-compose.yml` mounts `../models` and the AlphaFold cache read-only, so the
image stays small and weights can be regenerated or replaced without rebuilding.
For a public host, put nginx (or Caddy) in front of port 8765.

## Option C — Hugging Face Spaces (Docker SDK)

Create a Space with the Docker SDK, push this repository, set the root of the
Docker build to `deploy/` (or move the Dockerfile to the repository root), and
upload `models/` as a dataset mount or let the Space regenerate them at build
time. Spaces exposes port 7860 by default; change `EXPOSE` and the `--port`
argument accordingly.

## Health check and smoke test

```bash
curl -s http://HOST:8765/health | jq .
curl -s -X POST http://HOST:8765/predict \
  -H 'Content-Type: application/json' \
  -d '{"fasta":">P25101\nMEE...","threshold":0.9,"structure":true}' | jq '.results[0].sites[0:3]'
```

## Notes

- The service has no authentication or rate limiting; expose it behind a proxy
  that adds both if it is public.
- CPU-only is fine: scoring a 500-residue protein takes a few seconds; the model
  loads once and is reused.
- `--no-fetch` keeps the structure path offline-only, which is what a public
  deployment should use (no outbound calls from the server).
