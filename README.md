# Pareto optimization of masked superstrings

Supplementary repository for paper
*Pareto optimization of masked superstrings improves compression of pan-genome k-mer sets*. [**TODO** add bioRxiv link when possible]

This repository contains the implementation of Pareto optimization of masked superstrings for k-mer sets
and computation of lower bound for the number of runs of ones in the mask (or the number of matchtigs).

It also contains Snakemake pipelines reproducing the results from the paper,
and uses conda to manage software versions and most of the dependencies.

In case of any questions, feel free to ask by an [email](mailto:janci@kam.mff.cuni.cz).

## Howto

### Prepare

To install dependencies from conda and build localy provided software, use:
```bash
    ./prepare.sh
```

Next time, you only need to activate the conda environment with:
```bash
    conda activate ms-pareto-optimization
```

### Run

To run an experiment, move to the corresponding directory (`ex1-...`, `ex2-...`, `ex3-...`) and run:
```bash
    snakemake -j <number_of_threads> <any_optional_parameters>
```

Experiments produce `.tsv` files with results and plots
(**TODO** link new R plotting scripts to work from pipeline instead of manually).

### Modify

To tweak an experimental setup, modify the `Snakefile` in the corresponding directory.
Constants defined on top of Snakefiles define which datasets,
values of k, run penalties, and computation or compression methods are used.

## Project structure

Locally installed programs and snakefiles are stored in `tools` directory.

Other directory names can be modified in `tools/project_structure.smk` if needed.
Default names are:
- Datasets are downloaded into `data` directory.
- Computed superstring representations are stored in `computed` directory
  (it is recommended to turn off file search indexing for this directory).
- Results are stored in respective numbered directories (`ex1-...`, `ex2-...`, ...).

### Datasets

If you want to use the pipeline with custom datasets, you can modify `datasets.txt`.
Add the dataset name and url, the pipeline downloads datasets automatically
and can handle `xz`-compressed and uncompressed FASTA files.

In case you need to support other compression formats,
the simplest way is to modify the `download_data` rule in `tools/download.smk`.

## Implementation

Details about the implementation are provided in separate [README](./kmercamel-pareto/README.md).

The work was implemented inside a fork of KmerCamel🐫 and then copied over to another fork,
which results in the weird commit history with most of the relevant changes in the first commit made by @Jajopi
([d86ccf1a6b8007a6cb8680bb46f47a91ef058beb](https://github.com/Jajopi/ms-pareto-optimization-supplement/commit/d86ccf1a6b8007a6cb8680bb46f47a91ef058beb)).
The main concept (Pareto optimization) was internally called Joint optimization for a long time.
