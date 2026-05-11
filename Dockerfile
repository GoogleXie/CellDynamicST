# ============================================================================
# CellDynamicST Docker Image
# Full analysis environment with RStudio Server
# ============================================================================
FROM rocker/rstudio:4.4.0

LABEL maintainer="CellDynamicST Team"
LABEL description="CellDynamicST: Spatial Transcriptomic Analysis Pipeline"
LABEL version="0.1.0"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgeos-dev \
    libgdal-dev \
    libproj-dev \
    libglpk-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and Bioconductor dependencies
RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org')" \
    && R -e "BiocManager::install(c( \
        'org.Mm.eg.db', \
        'org.Hs.eg.db', \
        'GO.db', \
        'AnnotationDbi', \
        'clusterProfiler', \
        'enrichplot' \
    ), ask = FALSE, update = FALSE)"

# Install CRAN dependencies
RUN R -e "install.packages(c( \
    'Seurat', \
    'WGCNA', \
    'randomForest', \
    'mclust', \
    'igraph', \
    'pheatmap', \
    'viridis', \
    'ggrepel', \
    'patchwork', \
    'yaml', \
    'cli', \
    'checkmate', \
    'rlang', \
    'tidyr', \
    'dplyr', \
    'tibble', \
    'Matrix', \
    'future', \
    'future.apply', \
    'testthat', \
    'knitr', \
    'rmarkdown', \
    'quarto', \
    'remotes', \
    'renv', \
    'pkgdown', \
    'devtools' \
), repos = 'https://cran.r-project.org')"

# Copy package source
COPY . /home/rstudio/CellDynamicST
RUN chown -R rstudio:rstudio /home/rstudio/CellDynamicST

# Install the package
RUN R -e "devtools::install('/home/rstudio/CellDynamicST', dependencies = FALSE)"

# Copy notebooks to accessible location
RUN cp -r /home/rstudio/CellDynamicST/inst/notebooks /home/rstudio/notebooks \
    && cp /home/rstudio/CellDynamicST/inst/templates/experiment_config.yaml \
       /home/rstudio/notebooks/ \
    && chown -R rstudio:rstudio /home/rstudio/notebooks

# Create data directory
RUN mkdir -p /home/rstudio/data /home/rstudio/output \
    && chown -R rstudio:rstudio /home/rstudio/data /home/rstudio/output

# Expose RStudio Server port
EXPOSE 8787

# Default command: start RStudio Server
CMD ["/init"]
