# Snakemake workflow: Dada2

[![Snakemake](https://img.shields.io/badge/snakemake-≥5-brightgreen.svg)](https://snakemake.bitbucket.io)
[![dada2](https://img.shields.io/badge/dada2-v1.18-brightgreen.svg)](https://benjjneb.github.io/dada2/index.html)
<!-- [![Build Status](https://travis-ci.org/snakemake-workflows/amplicon-seq-dada2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/amplicon-seq-dada2) -->


This workflow is an implementation of the popular DADA2 tool. I followed the steps in the [Tutorial](https://benjjneb.github.io/dada2/tutorial.html). I utilized [Kraken2](https://ccb.jhu.edu/software/kraken2/) for ASV sequence classification instead of IDTaxa.

![dada2](https://benjjneb.github.io/dada2/images/DADA2_Logo_Text_1_14_640px.png)

The pipeline was inspired by the [Silas Kieser (@silask)'s dada2 snakemake pipeline](https://github.com/SilasK/16S-dada2).

## Authors

*  Rene Welch (@ReneWelch)

## Usage

### Install workflow


#### Install snakemake pipeline

Assuming [conda](https://docs.conda.io/en/latest/) is already installed, and a copy of this repository has been downloaded. Then, the pipeline can be installed by:

```sh
conda env create -n {env_name} --file dependencies.yml
```

#### Install R packages

This pipeline strongly depends on R package, therefore we utilized [renv](https://rstudio.github.io/renv/index.html) to set the right R package versions.

After `renv` is installed, utilize:

```R
renv::restore()
```

to install all the R libraries used in the pipeline

### Run workflow

* `snakemake -j{cores} all_profile` to plot the quality profiles
* `snakemake -j{cores}` to run the dada2 pipeline
* `snakemake -j{cores} all_taxonomy_kraken` to labels the ASV sequences with [kraken2](https://ccb.jhu.edu/software/kraken2/). Databases need to be downloaded from https://benlangmead.github.io/aws-indexes/k2
* `snakemake -j{cores} all_phyloseq` to perform the generate a `phyloseq` object, align sequences and compute the alpha / beta diversities. This step is not stable yet, but alternatively a good ending point could be `snakemake -j{cores} init_phyloseq`

* `snakemake -j{cores} clean` removes everything
* `snakemake -j{cores} clean_phyloseq` removes only the phyloseq object and diversities

## Cite

#### dada2

Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869

#### Kraken2:

Wood, D., Lu, J., Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biology 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0

#### phyloseq

McMurdie, P., Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS One 8, 4 (2013). https://doi.org/10.1371/journal.pone.0061217

