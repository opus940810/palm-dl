FROM python:3.12-slim

WORKDIR /app

# CPU-only scientific stack; the model binaries and AlphaFold cache are mounted
# at runtime (see docker-compose.yml) because they are large and regenerable.
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY src ./src
COPY server ./server
COPY data/processed ./data/processed

ENV PALMDL_THREADS=4
ENV PYTHONUNBUFFERED=1

EXPOSE 8765

# Default: sequence-only service. For the structure-aware path add
#   --structure --no-fetch  (requires a mounted results/structure/af cache)
CMD ["python", "server/app.py", "--host", "0.0.0.0", "--port", "8765"]
