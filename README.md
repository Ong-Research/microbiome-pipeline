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
**You may not need to do this step.**

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

## Set up data and metadata

1. **Set up a data/ directory with one subdirectory per batch.** The pipeline will look in the `data/` directory to find subdirectories that contain `fastq.gz` files. You can specify different filtering and trimming parameters per batch, and `dada2` will learn error rates per batch. So, you may wish to separate your files by sequencing batch and/or by sample type.
   
   ```sh
    mkdir data/
    mkdir data/batch01 [data/batch02 ... ]
   ```

   Next, populate the subdirectories with your `fastq.gz` (or `fastq`) files, making sure to include both R1 and R2. If your data are already located elsewhere in your file system, you can use symbolic links (symlinks) to avoid duplicating data:

   ```sh
    fqdir=/path/to/existing_fastq_files/
    batchdir=data/batch01
    ln -s ${fqdir}/*.fastq.gz ${batchdir}    
   ```

2. **Generate sample table.** The input for the pipeline are the sequencing files separated by batch in the `data/` directory. Using:

    ```sh
    Rscript ./prepare_sample_table.R
    ```

    will generate the `samples.tsv` file that contains 4 columns separated by tabs (the file will not actually contain |) 

    ```txt
    | batch | key | end1 | end2 |
    ```

    
3. **Set up sample metadata file.** This is where to link sample IDs to variables of interest to your study (condition, body site, time point, etc). **As currently implemented**, this is only used to augment the final TreeSummarizedExperiment file, which you will use for your own downstream analysis. You need to have a column `key` that matches the sample table. For example:


    ```txt
    | key | subject_id | body_site | time_point | ... 
    ```

    The default path for this file is `data/meta.tsv`. If you have it named otherwise, edit the variable `metadata` in `config/config.yaml` to point to your file.

4. **Set up negative control mapping file.** If you have negative control samples, you should make a table linking each study sample to its negative control kit(s). **As currently implemented**, this needs to be a `qs` file containing a tibble of this format:

    ```r
    # A tibble: 3 x 3
    batch        key       kits     
    <chr>        <chr>     <list>   
    1 batch01 sample_1  <chr [3]>
    2 batch01 sample_2  <chr [1]>
    3 batch02 sample_3  <chr [1]>
    ```

    Each element of `key` is a study sample ID matching the sample table.
    Each element of `kits` is a character list of sample keys for negative control samples,
    eg `c("sample_141", "sample_142", "sample_143")`.


    The default path for this file is `data/negcontrols.qs`. If you have it named otherwise, edit the variable `negcontroltable` in `config/config.yaml` to point to your file.
    
## Set up configuration file

Now you will edit `config.yaml` to set up running parameters. At the minimum, you will need to make these edits:

1. **Confirm general parameters and file paths.** Check that these variables are set to your satisfaction:

    ```sh
    # general parameters for run
    asv_prefix: "asv" # this will be the prefix for the ASV IDs
    sample_table: "samples.tsv" # this file is generated by prepare_sample_table.R
    negcontroltable: "data/negcontrols.qs" # manually generated
    metadata: "data/meta.tsv" # manually generated
    threads: 8 # max number of threads for a multithreaded snakemake rule, eg fastqc in quality_control
    ```

2. **Run sequence QC.** If you have not already looked at your data, the first thing you will want to do is run fastQC to generate quality profiles. 

    First, run `snakemake -j{cores} sequence_qc` to run fastQC, plot quality profiles, and build a `multiqc` report. You can find the individual plots in `workflow/report/quality_profiles` and the full report in `output/quality_control/multiqc`.

    Next, inspect those profiles. You may want to use them to decide where you want to trim R1 and R2 and how to set minimum overlap. 

3. **Edit parameters for each batch.** Check the `config.yaml` file for the default filter_and_trim, merge_pairs, and qc_parameters specifications for each batch. Keep in mind the length of your 16S region (eg, the V4 region is 252 bp long) and whether you will have enough overlap after trimming. Look up other `dada2` resources for more advice on how to choose these parameters.

    Here is an example for a 2x150 NextSeq run targeting the V4 region. This set of parameters keeps only full-length reads and requires an overlap of 25 base pairs. Because the V4 region is only ~252 base pairs long, we suspect that ASVs that are outside of that range may be bimeras. 

    ```txt
    filter_and_trim:
        batch01:
            truncQ: 0
            trimLeft: 0
            trimRight: 0
            maxLen: Inf
            minLen: 150
            maxN: 0
            minQ: 0
            maxEE: [2, 2]

    ...

    merge_pairs:
        batch01:
            minOverlap: 25
            maxMismatch: 0

    ...

    qc_parameters:
        negctrl_prop: 0.5
        max_length_variation: 5
        low_abundance_perc: 0.0001 
    ```


    The defaults provided in the pipeline are from a 2x300 run targeting V3-V4. Here we are truncating the ends of both reads because we noticed a significant loss of quality. 
    Because we have longer reads, we can require more overlap.

    ```txt
    filter_and_trim:
        dust_dec2018:
            truncQ: 12
            truncLen: [280, 250]
            trimLeft: 0
            trimRight: 0
            maxLen: Inf
            minLen: 100
            maxN: 0
            minQ: 0
            maxEE: [2, 2]
    ....

    merge_pairs:
        dust_dec2018:
        minOverlap: 50
        maxMismatch: 0

    ...

    qc_parameters:
        negctrl_prop: 0.5
        max_length_variation: 50
        low_abundance_perc: 0.0001 
    ```

## Download references for sequence classification

Databases need to be downloaded from <https://benlangmead.github.io/aws-indexes/k2>.
These are the four used by default in the config file. 

Here is an example of doing this using `wget`:

```
# download references
wdir="https://genome-idx.s3.amazonaws.com/kraken"
refdir="data/kraken_dbs"
mkdir -p $refdir

greengenes="16S_Greengenes13.5_20200326"
rdp="16S_RDP11.5_20200326"
silva="16S_Silva138_20200326"
minikraken="minikraken2_v2_8GB_201904"

for db in $greengenes $rdp $silva $minikraken;
do
    remote=${wdir}/${db}.tgz
    local=${refdir}/${db}.tgz
    wget $remote -O $local
    tar -xzf ${local} -C ${refdir}
done
```

## Workflow commands

Replace `{cores}` with the maximum number of cores that you want
to be used at one time.

The rules `taxonomy` and `phylotree` will take a little longer the first
time they are run because they install conda environments.

* `snakemake -j{cores}` runs everything
* `snakemake -j{cores} sequence_qc` plots quality profiles, and build a `multiqc` report
* `snakemake -j{cores} [--nt] dada2` computes the ASV matrix from the different batches (use --nt option to keep intermediate files for troubleshooting)
* `snakemake -j{cores} taxonomy --use-conda` to labels the ASV sequences with [kraken2](https://ccb.jhu.edu/software/kraken2/). Databases need to be downloaded in advance from <https://benlangmead.github.io/aws-indexes/k2>
* `snakemake -j{cores} phylotree --use-conda` computes the phylogenetic tree using `qiime2`'s FastTree
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