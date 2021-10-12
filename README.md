# Snakemake workflow: Dada2

[![Snakemake](https://img.shields.io/badge/snakemake-≥5-brightgreen.svg)](https://snakemake.bitbucket.io)
[![dada2](https://img.shields.io/badge/dada2-v1.18-brightgreen.svg)](https://benjjneb.github.io/dada2/index.html)
<!-- [![Build Status](https://travis-ci.org/snakemake-workflows/amplicon-seq-dada2.svg?branch=master)](https://travis-ci.org/snakemake-workflows/amplicon-seq-dada2) -->

This workflow is an implementation of the popular DADA2 tool. I followed the steps in the [Tutorial](https://benjjneb.github.io/dada2/tutorial.html). I utilized [Kraken2](https://ccb.jhu.edu/software/kraken2/) for ASV sequence classification instead of IDTaxa.

![dada2](https://benjjneb.github.io/dada2/images/DADA2_Logo_Text_1_14_640px.png)

The pipeline was inspired by the [Silas Kieser (@silask)'s dada2 snakemake pipeline](https://github.com/SilasK/16S-dada2).

## Authors

* Rene Welch (@ReneWelch)

## Install workflow

There are two steps to install the pipeline:

1. **Install snakemake and dependencies:** Assuming [conda](https://docs.conda.io/en/latest/) is already installed, and a copy of this repository has been downloaded. Then, the pipeline can be installed by:

    ```sh
    conda env create -n {env_name} --file dependencies.yml
    conda env create -n microbiome --file dependencies.yml
    ```

2. **Install R packages:** To install the R packages used by this pipeline use:

    ```sh
    R CMD BATCH --vanilla ./install_r_packages.R
    ```

### Troubleshooting conda and environment variables

If you have other versions of R and R user libraries elsewhere, 
you might encounter some problems with environment variables and conda. 
You may need to provide a local `.bashrc` file to 
place the conda path at the beginning of your `$PATH` environment variable.

    ```sh
    source /etc/bash_completion.d/git
    __conda_setup="$('/path/to/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else
        if [ -f "/path/to/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/path/to/miniconda3/etc/profile.d/conda.sh"
        else
            export PATH="/path/to/miniconda3/bin:$PATH"
        fi
    fi
    unset __conda_setup

    alias R='R --vanilla'
    export PATH="/path/to/miniconda3/envs/microbiome/bin:$PATH"
    ```

## Set up data directories and metadata files

1. **Set up a data/ directory with one subdirectory per batch.** The pipeline will look in the `data/` directory to find subdirectories that contain `fastq.gz` files. You can specify different filtering and trimming parameters per batch, and `dada2` will learn error rates per batch. So, you may wish to separate your files by sequencing batch and/or by sample type.
   
   ```sh
    mkdir data/
    mkdir data/batch01 \[data/batch02 ... \]
   ```

   Next, populate the subdirectories with your `fastq.gz` (or `fastq`) files, making sure to include both R1 and R2. If your data are already located elsewhere in your file system, you can use symbolic links (symlinks) to avoid duplicating data:

   ```sh
    fqdir=/path/to/existing_fastq_files/
    batchdir=/data/batch01
    ln -s ${fqdir}/*.fastq.gz ${batchdir}    
   ```

2. **Generate sample table.** The input for the pipeline are the sequencing files separated by batch in the `data/` directory. Using:

    ```sh
    Rscript ./prepare_sample_table.R
    ```

will generate the `samples.tsv` file that contains 4 columns (the file will not actually contain |)

    ```txt
    | batch | key | end1 | end2 |
    ```

separated by a tab space. 

3. **Set up sample metadata file.**

4. **Set up negative control mapping file.** If you have negative control samples, you should make a table linking each study sample to its negative control kit(s). **At this point**, this needs to be a `qs` file containing a tibble of this format:

    ```r
        # A tibble: ... x 3
    batch        key       kits     
    <chr>        <chr>     <list>   
    1 batch01 sample_1  <chr [3]>
    2 batch01 sample_2  <chr [1]>
    3 batch02 sample_3  <chr [1]>
    ```

Each element of `key` is a study sample ID matching the sample table.
Each element of `kits` is a character list of sample keys for negative control samples,
eg `c("sample_141", "sample_142", "sample_143")`.
    
## Run workflow

Then, we can use different commands in the pipeline, for example considering 16 threads:

* `snakemake -j{cores}` runs everything
* `snakemake -j{cores} sequence_qc` plots quality profiles, and build a `multiqc` report
* `snakemake -j{cores} dada2` computes the ASV matrix from the different batches
* `snakemake -j{cores} all_taxonomy_kraken` to labels the ASV sequences with [kraken2](https://ccb.jhu.edu/software/kraken2/). Databases need to be downloaded from <https://benlangmead.github.io/aws-indexes/k2>
* `snakemake -j{cores} phylotree` computes the phylogenetic tree using `qiime2`'s FastTree
* `snakemake -j{cores} mia` prepare the `TreeSummarizedExperiment` containing all the data generated

![microbiome_pipeline](microbiome.png)

## Cite

### dada2

Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). <https://doi.org/10.1038/nmeth.3869>

### Kraken2

Wood, D., Lu, J., Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biology 20, 257 (2019). <https://doi.org/10.1186/s13059-019-1891-0>

### phyloseq

McMurdie, P., Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS One 8, 4 (2013). <https://doi.org/10.1371/journal.pone.0061217>

### qiime2

Bolyen, Evan, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian Abnet, Gabriel A. Al-Ghalith, Harriet Alexander, et al. 2018. “QIIME 2: Reproducible, Interactive, Scalable, and Extensible Microbiome Data Science.” e27295v2. PeerJ Preprints. <https://doi.org/10.7287/peerj.preprints.27295v2>.

### FastTree

Price, Morgan N., Paramvir S. Dehal, and Adam P. Arkin. 2010. “FastTree 2--Approximately Maximum-Likelihood Trees for Large Alignments.” PloS One 5 (3): e9490.

### TreeSummarizedExperiment

Huang, Ruizhu, Charlotte Soneson, Felix G. M. Ernst, Kevin C. Rue-Albrecht, Guangchuang Yu, Stephanie C. Hicks, and Mark D. Robinson. 2020. “TreeSummarizedExperiment: A S4 Class for Data with Hierarchical Structure.” F1000Research 9 (October): 1246.