############################################################
## Dockerfile for StrainQuant
############################################################

FROM ubuntu:22.04
LABEL maintainer="Adam Dinan <adam1989ie@gmail.com>"

ARG DEBIAN_FRONTEND=noninteractive

ENV PYTHON_VERSION=3.8.5

ADD https://raw.githubusercontent.com/adamd3/StrainQuant/main/requirements.txt .

COPY requirements.txt /tmp
WORKDIR /tmp


RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates \
    build-essential \
    curl libssl-dev libcurl4-openssl-dev libxml2-dev \
    python3-numpy python3-pip gawk pigz r-base-dev fastqc \
    trim-galore kallisto && \
    apt-get clean autoclean

RUN pip install --no-cache-dir -r requirements.txt

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.rstudio.com'; \
    options(repos = r);" > ~/.Rprofile

RUN R -e 'install.packages(c(  \
    "optparse", "RColorBrewer", "reshape2", "tidyverse",   \
    "readr", "tibble", "purrr", "tidyr", "stringr", "ggplot2", \
    "scales", "pheatmap", "matrixstats", "umap", "BiocManager"))'

RUN R -e 'BiocManager::install(c("edgeR", "DESeq2"))'
