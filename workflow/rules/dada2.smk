# An example collection of Snakemake rules imported in the main Snakefile.
    
rule filter_and_trim:
  input:
    R1= SampleTable.R1.values,
    R2= SampleTable.R2.values
  output:
    R1 = expand("data/filtered/{sample}_R1.fastq.gz",sample=SAMPLES),
    R2 = expand("data/filtered/{sample}_R2.fastq.gz",sample=SAMPLES),
    nreads = temp("data/stats/Nreads_filtered.txt")
  params:
    samples = SAMPLES
  threads:
    config["threads"]
  # conda:
  #   "../envs/dada2.yaml"
  log:
    "logs/dada2/filter.txt"
  script:
    "../scripts/dada2/filter_and_trim.R"


rule learn_error_rates:
  input:
    R1 = rules.filter_and_trim.output.R1,
    R2 = rules.filter_and_trim.output.R2
  output:
    errR1 = "data/model/error_rates_R1.qs",
    errR2 = "data/model/error_rates_R2.qs",
    plotErr1 = "figures/model/error_rates_R1.png",
    plotErr2 = "figures/model/error_rates_R2.png"
  threads:
    config["threads"]
  # conda:
  #     "../envs/dada2.yaml"
  log:
    "logs/dada2/learnErrorRates.txt"
  script:
    "../scripts/dada2/learn_error_rates.R"

rule dereplicate_single:
  input:
    R1 = "data/filtered/{sample}_R1.fastq.gz",
    R2 = "data/filtered/{sample}_R2.fastq.gz",
    errR1 = rules.learn_error_rates.output.errR1,
    errR2 = rules.learn_error_rates.output.errR2
  output:
    merge = temp("data/mergers/{sample}_asv.qs")
  log:
    "logs/dada2/{sample}_merge.txt"
  script:
    "../scripts/dada2/dereplicate_one_sample_pair.R"

rule dereplicate:
  input:
    derep = expand("data/mergers/{sample}_asv.qs", sample = SAMPLES)
  output:
    seqtab = temp("data/asv/seqtab_with_chimeras.qs"),
    nreads = temp("data/stats/Nreads_dereplicated.txt")
  params:
    samples = SAMPLES
  threads:
    config["threads"]
  log:
    "logs/dada2/dereplicate_samples.txt"
  script:
    "../scripts/dada2/gather_derep_seqtab.R"


rule remove_chimeras:
  input:
    seqtab = rules.dereplicate.output.seqtab
  output:
    asvs = temp("data/asv/seqtab_nochimeras.qs"),
    nreads = "data/stats/Nreads_chimera_removed.txt"
  threads:
    config["threads"]
  # conda:
  #       "../envs/dada2.yaml"
  log:
    "logs/dada2/remove_chimeras.txt"
  script:
    "../scripts/dada2/remove_chimeras.R"

rule filter_asvs:
  input:
    seqtab = rules.remove_chimeras.output.asvs,
    negcontrol = config["negcontroltable"]
  output:
    plot_seqlength = "figures/qc/nasvs_by_seqlength.png",
    plot_seqabundance = "figures/qc/nasvs_by_seqabundance.png",
    seqtab_filt = "data/asv/seqtab_nochimeras_qc.qs"
  log:
    "logs/dada2/filter_qc.txt"
  script:
    "../scripts/dada2/filter_asvs.R"

rule stats:
  input:
    nreads_filtered = "data/stats/Nreads_filtered.txt",
    nreads_dereplicated = "data/stats/Nreads_dereplicated.txt",
    nreads_chim_removed = "data/stats/Nreads_chimera_removed.txt"
  output:
    nreads = "data/stats/Nreads_dada2.txt",
    fig_step = "figures/qc/dada2steps_vs_abundance.png",
    fig_step_rel = "figures/qc/dada2steps_vs_relabundance.png",
  log:
    "logs/dada2/summarize_stats.txt"
  script:
    "../scripts/dada2/summarize_nreads.R"

# rule filterLength:
#     input:
#         seqtab= rules.removeChimeras.output.rds
#     output:
#         plot_seqlength= "figures/Lengths/Sequence_Length_distribution.pdf",
#         plot_seqabundance= "figures/Lengths/Sequence_Length_distribution_abundance.pdf",
#         rds= "output/seqtab.rds",
#         tsv=  "output/seqtab.tsv",
#     threads:
#         1
#     conda:
#         "../envs/dada2.yaml"
#     log:
#         "logs/dada2/filterLength.txt"
#     script:
#         "../scripts/dada2/filterLength.R"



# rule IDtaxa:
#     input:
#         seqtab= "output/seqtab.rds",
#         ref= lambda wc: config['idtaxa_dbs'][wc.ref]
#     output:
#         taxonomy= "taxonomy/{ref}.tsv",
#     threads:
#         config["threads"]
#     log:
#         "logs/dada2/IDtaxa_{ref}.txt"
#     script:
#         "../scripts/dada2/IDtaxa.R"

# localrules: get_ggtaxonomy

# rule get_ggtaxonomy:
#     input:
#         rules.IDtaxa.output
#     output:
#         taxonomy= "taxonomy/{ref}_gg.tsv",
#     threads:
#         1
#     log:
#         "logs/dada2/get_ggtaxonomy_{ref}.txt"
#     run:

#         import pandas as pd

#         tax = pd.read_csv(input[0],sep='\t',index_col=0)

#         out= tax.dropna(how='all',axis=1).dropna(how='all',axis=0)

#         out= out.apply(lambda col: col.name[0]+'_'+col.dropna())
#         out= out.apply(lambda row: ';'.join(row.dropna()),axis=1)
#         out.name='Taxonomy'


#         out.to_csv(output[0],sep='\t',header=True)







# rule get_rep_seq:
#     input:
#         "output/seqtab.tsv",
#     output:
#         "taxonomy/rep_seq.fasta"
#     run:
#         with open(input[0]) as infile:
#             seqs = infile.readline().strip().replace('"','').split()
#         with open(output[0],'w') as outfile:
#             for s in seqs:
#                 outfile.write(f">{s}\n{s}\n")
