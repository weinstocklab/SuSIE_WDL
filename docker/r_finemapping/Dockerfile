FROM rocker/r-ver:4.5.3

# System libraries required by R packages:
#   Rfast  -> BLAS/LAPACK + Fortran
#   ggplot2 ggsave (PDF/PNG) -> Cairo, fontconfig, freetype
RUN apt-get update && apt-get install -y --no-install-recommends \
        liblapack-dev \
        libopenblas-dev \
        gfortran \
        libcairo2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Pre-install all R packages used by run_susie and plot_locus
RUN Rscript -e "\
    install.packages( \
        c('susieR', 'data.table', 'Rfast', 'ggplot2', 'ggrepel', 'cowplot'), \
        repos = 'https://cloud.r-project.org', \
        quiet = TRUE \
    )"

# Bundle gencode annotation data (from weinstocklab/gwasplot) so
# plot_locus does not need a network download at runtime.
RUN curl -fsSL \
    https://raw.githubusercontent.com/weinstocklab/gwasplot/main/data/gencode.rda \
    -o /usr/local/lib/R/gencode.rda
