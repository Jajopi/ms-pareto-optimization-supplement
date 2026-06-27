#!/bin/bash

set -euo pipefail

# Almost all dependencies are installed through conda, except for the ones listed below in this file:
conda config --add channels conda-forge
conda config --add channels bioconda
conda install GeCo3 ggcat bioconda::jellyfish kmercamel snakemake xz

# Those dependencies are not distributed through conda, so they are included in the repository:

# Our-work - Pareto optimization, matchtig count lowerbound, sorted integers
make -C kmercamel-pareto

# Elias-Fano encoding for ex3
{
  cd tools/EliasFano/2i_bench; cmake .; cd ../../../
  make -C tools/EliasFano/2i_bench
  make -C tools/EliasFano
} || true # This thing is old and weird and I was only able to get it compiled to 25%, but that's apparently enough
