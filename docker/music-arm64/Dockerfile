FROM rocker/r-ver:4.5.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('BiocManager','remotes','peakRAM'))"
RUN R -q -e "options(repos = BiocManager::repositories()); BiocManager::install(c('Biobase','IRanges','S4Vectors','S4Arrays','DelayedArray','GenomeInfoDbData','GenomeInfoDb','GenomicRanges','SummarizedExperiment','SingleCellExperiment'), ask = FALSE, update = FALSE); for (p in c('Biobase','IRanges','S4Vectors','S4Arrays','DelayedArray','GenomeInfoDbData','GenomeInfoDb','GenomicRanges','SummarizedExperiment','SingleCellExperiment')) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf('%s install failed', p))"
RUN R -q -e "options(repos = BiocManager::repositories()); BiocManager::install('TOAST', ask = FALSE, update = FALSE); if (!requireNamespace('TOAST', quietly = TRUE)) stop('TOAST install failed')"
RUN R -q -e "options(repos = BiocManager::repositories()); BiocManager::install('MuSiC', ask = FALSE, update = FALSE); if (!requireNamespace('MuSiC', quietly = TRUE)) remotes::install_github('xuranw/MuSiC', upgrade = 'never'); if (!requireNamespace('MuSiC', quietly = TRUE)) stop('MuSiC install failed')"

WORKDIR /workspace
