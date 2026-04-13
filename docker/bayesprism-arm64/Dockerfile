FROM rocker/r-ver:4.5.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libglpk40 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('remotes','peakRAM','snowfall','gplots','BiocManager'))"
RUN R -q -e "options(repos = BiocManager::repositories()); BiocManager::install(c('Biobase','XVector','SparseArray','DelayedArray','SummarizedExperiment','SingleCellExperiment','BiocParallel','bluster'), ask = FALSE, update = FALSE); if (!requireNamespace('Biobase', quietly = TRUE)) stop('Biobase install failed'); if (!requireNamespace('XVector', quietly = TRUE)) stop('XVector install failed'); if (!requireNamespace('bluster', quietly = TRUE)) stop('bluster install failed')"
RUN R -q -e "options(repos = BiocManager::repositories()); BiocManager::install('scran', ask = FALSE, update = FALSE); if (!requireNamespace('scran', quietly = TRUE)) stop('scran install failed')"
RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages('NMF'); if (!requireNamespace('NMF', quietly = TRUE)) stop('NMF install failed')"
RUN git clone --depth 1 https://github.com/Danko-Lab/BayesPrism.git /tmp/BayesPrism \
    && R -q -e "remotes::install_local('/tmp/BayesPrism/BayesPrism', dependencies = TRUE, upgrade = 'never')" \
    && rm -rf /tmp/BayesPrism

WORKDIR /workspace
