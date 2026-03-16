#!/bin/bash
set -euo pipefail
snakemake --rerun-incomplete $*
