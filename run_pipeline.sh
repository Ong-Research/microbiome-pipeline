# commands run in order
# DC 2021-10-11

# don't forget to activate!!
#source .bashrc
#conda activate microbiome

# quality control
#snakemake -j8 sequence_qc

# dada2 to generate ASV table
#snakemake -j10 --nt dada2
# --nt will keep  temp files if you want to troubleshoot
 
# download references
#wdir="https://genome-idx.s3.amazonaws.com/kraken"
#refdir="data/kraken_dbs"
#mkdir -p $refdir

#greengenes="16S_Greengenes13.5_20200326"
#rdp="16S_RDP11.5_20200326"
#silva="16S_Silva138_20200326"
#minikraken="minikraken2_v2_8GB_201904"
#for db in $greengenes $rdp $silva $minikraken;
#do
#    remote=${wdir}/${db}.tgz
#    local=${refdir}/${db}.tgz
#    echo wget $remote -O $local
#done

#snakemake -j8 taxonomy --use-cond # to labels the ASV sequences with [kraken2](https://ccb.jhu.edu/software/kraken2/). 
# Databases need to be downloaded in advance from <https://benlangmead.github.io/aws-indexes/k2>
#snakemake -j8 phylotree --use-conda # computes the phylogenetic tree using `qiime2`'s FastTree
#snakemake -j8 mia # prepare the `TreeSummarizedExperiment` containing all the data generated
