#!/bin/bash

set -euo pipefail

ENV_NAME="ms-pareto-optimization"

# Almost all dependencies are installed through conda, except for the ones listed below in this file
if conda env list | grep "$ENV_NAME" > /dev/null 2>&1; then
    echo "Environment "$ENV_NAME" already present."
    conda activate "$ENV_NAME"
else
    conda env create --name "$ENV_NAME"
    conda activate "$ENV_NAME"
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install snakemake kmercamel xz GeCo3 matchtigs bioconda::jellyfish
    pip install matplotlib
fi

# Pareto optimization, matchtig count lowerbound, sorted integers
make -C kmercamel-pareto

# Elias-Fano encoding
make -C tools/EliasFano/2i_bench
make -C tools/EliasFano

# Optimal-length matchtigs
make -C tools/blossom5
