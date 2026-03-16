# Howto

To install dependencies from conda and localy provided programs, use:
```bash
    ./prepare.sh
```

Next time, you only need to activate the conda environment with:
```bash
    conda activate ms-pareto-optimization
```

To run an experiment, move to the corresponding directory (`ex1-...`, `ex2-...`, `ex3-...`) and run:
```bash
    snakemake -j <number_of_threads> <any_optional_parameters>
```

# Project structure

Locally installed programs and snakefiles are stored in `tools` directory.

Other directory names can be modified in `tools/project_structure.smk` if needed.
Default names are:
- Datasets are downloaded into `data` directory.
- Computed superstring representations are stored in `computed` directory
  (it is recommended to turn off file search indexing for this directory).
- Results are stored in respective numbered directories (`ex1-...`, `ex2-...`, ...).

## Datasets

If you want to use the pipeline with custom datasets, you can modify `datasets.txt`.
Add the dataset name and url, the pipeline downloads datasets automatically
and can handle `xz`-compressed and uncompressed FASTA files.

In case you need to support other compression formats,
the simplest way is to modify the `download_data` rule in `tools/download.smk`.
