# commands run in order
# DC 2021-10-11

# don't forget to activate!!
source .bashrc
conda activate microbiome

# quality control
#snakemake -j8 sequence_qc

# dada2 to generate ASV table
snakemake -j10 dada2 