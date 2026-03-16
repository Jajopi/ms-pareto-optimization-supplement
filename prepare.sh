#!/bin/bash

set -euo pipefail

# Almost all dependencies are installed through conda, except for the ones listed below in this file
if conda env list | grep 'ms-pareto-optimization-supplement' > /dev/null 2>&1; then; else conda env create --file environment.yml; fi
conda activate ms-pareto-optimization-supplement

# Pareto optimization, matchtig count lowerbound
make -C kmercamel

# Elias-Fano encoding
make -C tools/EliasFano/2i_bench
make -C tools/EliasFano

# optimal-length matchtigs
make -C tools/blossom5
