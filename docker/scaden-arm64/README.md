# Scaden ARM64 Container

Build an ARM64-compatible Scaden image locally:

```bash
docker build -t rectangle-scaden:arm64 -f docker/scaden-arm64/Dockerfile docker/scaden-arm64
```

Optional explicit platform build:

```bash
docker build --platform linux/arm64/v8 -t rectangle-scaden:arm64 -f docker/scaden-arm64/Dockerfile docker/scaden-arm64
```

The benchmark script `notebooks/scaden.py` defaults to:
- `SCADEN_DOCKER_IMAGE=rectangle-scaden:arm64`
- `SCADEN_DOCKER_PLATFORM=linux/arm64/v8`

You can override either with environment variables.
