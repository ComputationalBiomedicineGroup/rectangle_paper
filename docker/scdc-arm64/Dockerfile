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
    libsqlite3-dev \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('BiocManager','remotes','peakRAM'))"
RUN R -q -e "options(repos = BiocManager::repositories()); BiocManager::install(c('Biobase','BiocGenerics','Biostrings','AnnotationDbi','KEGGREST'), ask = FALSE, update = FALSE); for (p in c('Biobase','BiocGenerics','Biostrings','AnnotationDbi','KEGGREST')) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf('%s install failed', p))"
RUN R -q -e "remotes::install_github('renozao/xbioc', upgrade = 'never'); if (!requireNamespace('xbioc', quietly = TRUE)) stop('xbioc install failed')"
RUN R -q -e "remotes::install_github('meichendong/SCDC', upgrade = 'never'); if (!requireNamespace('SCDC', quietly = TRUE)) stop('SCDC install failed')"

WORKDIR /workspace
