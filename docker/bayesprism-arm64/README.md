# BayesPrism ARM64 Container

Build an ARM64-compatible BayesPrism image locally:

```bash
docker build -t rectangle-bayesprism:arm64 -f docker/bayesprism-arm64/Dockerfile docker/bayesprism-arm64
```

Optional explicit platform build:

```bash
docker build --platform linux/arm64/v8 -t rectangle-bayesprism:arm64 -f docker/bayesprism-arm64/Dockerfile docker/bayesprism-arm64
```

Run the benchmark from repo root:

```bash
docker run --rm \
  --platform linux/arm64/v8 \
  -v "$PWD":/workspace \
  -w /workspace \
  rectangle-bayesprism:arm64 \
  Rscript notebooks/bayesprism.R
```

Optional custom inputs:

```bash
docker run --rm \
  --platform linux/arm64/v8 \
  -v "$PWD":/workspace \
  -w /workspace \
  -e HAO_SCE_PATH=/workspace/data/Hao_sce.rds \
  -e BULKS_PATH=/workspace/data/bulks_finotello_counts.csv \
  rectangle-bayesprism:arm64 \
  Rscript notebooks/bayesprism.R
```

`notebooks/bayesprism.R` writes `bayesprism_benchmark_results.csv` into `/workspace`.
