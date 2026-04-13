FROM rocker/r-ver:4.5.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('BiocManager','remotes','peakRAM'))"
RUN R -q -e "BiocManager::install('Biobase', ask = FALSE, update = FALSE)"
RUN R -q -e "remotes::install_github('cozygene/bisque', upgrade = 'never')"

WORKDIR /workspace
