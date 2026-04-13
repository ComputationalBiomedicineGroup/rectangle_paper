# Bisque ARM64 Container

Build an ARM64-compatible Bisque image locally:

```bash
docker build -t rectangle-bisque:arm64 -f docker/bisque-arm64/Dockerfile docker/bisque-arm64
```

Optional explicit platform build:

```bash
docker build --platform linux/arm64/v8 -t rectangle-bisque:arm64 -f docker/bisque-arm64/Dockerfile docker/bisque-arm64
```

Run the benchmark from repo root:

```bash
docker run --rm \
  --platform linux/arm64/v8 \
  -v "$PWD":/workspace \
  -w /workspace \
  rectangle-bisque:arm64 \
  Rscript notebooks/bisque.R
```

Optional custom inputs:

```bash
docker run --rm \
  --platform linux/arm64/v8 \
  -v "$PWD":/workspace \
  -w /workspace \
  -e HAO_SCE_PATH=/workspace/data/Hao_sce.rds \
  -e BULKS_PATH=/workspace/data/bulks_finotello_tpm.csv \
  rectangle-bisque:arm64 \
  Rscript notebooks/bisque.R
```

`notebooks/bisque.R` writes `bisque_benchmark_results.csv` into the working directory (`/workspace`).
