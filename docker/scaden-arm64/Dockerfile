FROM python:3.10-slim

ENV PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    TF_CPP_MIN_LOG_LEVEL=2

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade pip setuptools wheel

# Pin versions for reproducible benchmarking on ARM64.
RUN pip install \
    "setuptools<81" \
    "numpy<2" \
    "tensorflow==2.15.1" \
    "scaden==1.1.2"

WORKDIR /workspace
