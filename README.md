# Pareto optimization of masked superstrings

Supplementary repository for paper
[Pareto optimization of masked superstrings improves compression of pan-genome k-mer sets](https://www.biorxiv.org/content/10.64898/2026.03.18.712440).

This repository contains the implementation of *Pareto optimization of masked superstrings for k-mer sets*
and computation of lower bound for the number of runs of ones in the mask (or the number of matchtigs).

It also contains Snakemake pipelines reproducing the results from the paper,
and uses conda to manage software versions and most of the dependencies.

We plan to also add the exact scripts creating the figures in the paper; however, this is currently still a TODO.

In case of any questions, feel free to ask by an [email](mailto:janci@kam.mff.cuni.cz).

## Howto

### Prepare

To install dependencies from conda and build localy provided software, use:
```bash
    conda env create --name ms-pareto-optimization
    conda activate ms-pareto-optimization
    ./prepare.sh
```

Next time, you only need to activate the conda environment with:
```bash
    conda activate ms-pareto-optimization
```

### Run

To run an experiment, move to the corresponding directory (`ex1-...`, `ex2-...`, ...) and run:
```bash
    rm -rf results/<dataset for which you want to recompute the results>
    snakemake -j <number_of_threads_to_use> <any_optional_parameters_of_snakemake>
```

Experiments produce `.tsv` files with results (and, in the future, plots).
If you used `git clone` to obtain the repository, you can compare them to the original results using `git diff results`. 

Example snakemake configuration: `snakemake -j 10 --rerun-incomplete --resources mem_mb=26000`

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
- Results are stored in respective numbered directories (`ex1-...`, `ex2-...`, ...) in a `results` subdirectories.

### Datasets

If you want to use the pipeline with custom datasets, you can modify `datasets.txt`.
Add the dataset name and url, the pipeline downloads datasets automatically
and can handle `xz`-compressed and uncompressed FASTA files.

In case you need to support other compression formats,
you can modify the `download_data` rule in `tools/download.smk`,
but probably the simplest way is to just download and extract the file by yourself,
in which case the downloading part of the pipeline is not used at all.

## Implementation

Details about the implementation are linked in separate [README](./kmercamel-pareto/README.md).

The work was implemented inside a fork of KmerCamel🐫 and then copied over to another fork,
which resulted in a weird commit history cor a while.
Now it is settled like this:
- the first commit contains the history of original KmerCamel🐫, which has been forked.
- the relevant changes (Pareto optimization implementation) are in the second commit made by @Jajopi
([d86ccf1a6b8007a6cb8680bb46f47a91ef058beb](https://github.com/Jajopi/ms-pareto-optimization-supplement/commit/d86ccf1a6b8007a6cb8680bb46f47a91ef058beb)).

The main concept (Pareto optimization) was internally called Joint optimization for a long time.
All occurences in the files should be up to date; however, you may encounter the old name in the history.
