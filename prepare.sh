#!/bin/bash

set -euo pipefail

# Almost all dependencies are installed through conda, except for the ones listed below in this file:
conda config --add channels conda-forge
conda config --add channels bioconda
conda install GeCo3=1.0 ggcat=2.0.0 bioconda::jellyfish=2.2.10 kmercamel=2.2.0 snakemake=9.22.0 xz=5.8.1 bzip2=1.0.8

# Those dependencies are not distributed through conda, so they are included in the repository:

# Our-work - Pareto optimization, matchtig count lowerbound, sorted integers
make -C kmercamel-pareto

# Elias-Fano encoding for ex3
{
  cd tools/EliasFano/2i_bench; rm CMakeCache.txt; cmake . -DCMAKE_POLICY_VERSION_MINIMUM=3.5; cd ../../../
  make -C tools/EliasFano/2i_bench # This thing is old and weird and I was only able to get it compiled to 25%, but that's apparently enough
  make -C tools/EliasFano
} || true
