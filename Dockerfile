FROM condaforge/miniforge3:latest

RUN mamba install -y -c conda-forge -c bioconda \
        r-base=4.4.2 \
        r-data.table=1.17.8 \
        r-r.utils=2.13.0 \
        bcftools=1.23 \
        bedtools=2.31.1 \
        htslib=1.23 \
    && mamba clean -afy

LABEL maintainer="Mark" \
      description="DMRdetector ONT Methylation Analysis Pipeline container"
