rule filter_and_trim_sample:
  input:
    end1 = lambda wc: sample_dict[wc.sample]["end1"],
    end2 = lambda wc: sample_dict[wc.sample]["end2"]
  output:
    end1 = "output/dada2/filtered/{sample}_filtered_R1.fastq.gz",
    end2 = "output/dada2/filtered/{sample}_filtered_R2.fastq.gz",
    summary = "output/dada2/summary/{sample}_summary_filtered.tsv"
  params:
    config = "config/config.yaml",
    sample_name = "{sample}",
    batch = lambda wc: sample_dict[wc.sample]["batch"]
  log:
    "logs/dada2/01_filter_and_trim_{sample}.log"
  shell:
    """Rscript workflow/scripts/dada2/filter_and_trim.R \
      {output.end1} {output.end2} {output.summary} \
      {params.sample_name} --end1={input.end1} --end2={input.end2} \
      --batch={params.batch} --log={log} --config={params.config}"""

rule learn_error_rates_batch_end:
  input:
    filtered = lambda wc: expand("output/dada2/filtered/{sample}_filtered_{end}.fastq.gz", sample = get_batch_samples(sample_table, wc.batch), end = ["R1", "R2"])
  output:
    mat = "output/dada2/model/{batch}_error_rates_{end}.qs",
    plot = "workflow/report/model/{batch}_error_rates_{end}.png"
  threads: config["threads"] / 2
  params:
    # sample_names = lambda wc: get_batch_samples(sample_table, wc.batch),
    config = "config/config.yaml",
    batch = "{batch}"
  log:
    "logs/dada2/02_learn_error_rates_{batch}_{end}.log"
  shell:
    """Rscript workflow/scripts/dada2/learn_error_rates.R \
      --error_rates={output.mat} --plot_file={output.plot} \
      {input.filtered} --log={log} --batch={params.batch} \
      --config={params.config} --cores={threads}"""

def get_filtered_files(wc, end):
  """
  gets the name of the filtered files based on a wildcard
  """
  out = "output/dada2/filtered/" + wc.sample + "_filtered_" + end + ".fastq.gz"
  return out

def get_error_matrix(wc, sample_table, end):
  """
  gets the error matrix based on a wildcard
  """
  sample_dict = sample_table.to_dict("index")
  batch = sample_dict[wc.sample]["batch"]
  out = "output/dada2/model/" + batch + "_error_rates_" + end + ".qs"
  return out

rule dereplicate_sample:
  input:
    filt_end1 = lambda wc: get_filtered_files(wc, "R1"),
    filt_end2 = lambda wc: get_filtered_files(wc, "R2"),
    mat_end1 = lambda wc: get_error_matrix(wc, sample_table, "R1"),
    mat_end2 = lambda wc: get_error_matrix(wc, sample_table, "R2")
  output:
    merge = "output/dada2/merge/{batch}/{sample}_asv.qs"
  params:
    config = "config/config.yaml",
    sample_name = lambda wc: wc.sample,
    batch = lambda wc: sample_dict[wc.sample]["batch"]
  threads: 1
  log:
    "logs/dada2/03_{batch}_{sample}_merge.txt"
  shell:
    """Rscript workflow/scripts/dada2/dereplicate_one_sample_pair.R \
      {output.merge} {params.sample_name} \
      {input.filt_end1} {input.filt_end2} \
      --end1_err={input.mat_end1} --end2_err={input.mat_end2} \
      --log={log} --batch={params.batch} --config={params.config}"""

rule dereplicate_batch:
  input:
    derep = lambda wc: expand("output/dada2/merge/{batch}/{sample}_asv.qs", sample = get_batch_samples(sample_table, wc.batch), batch = wc.batch)
  output:
    asv = "output/dada2/asv_batch/{batch}_asv.qs",
    summary = "output/dada2/asv_batch/{batch}_summary.qs"
  params:
    config = "config/config.yaml",
    batch = lambda wc: wc.batch
  threads: 1
  log:
    "logs/dada2/04_{batch}_derep_batch.log"
  shell:
    """Rscript workflow/scripts/dada2/gather_derep_seqtab.R \
      {output.asv} {output.summary} {input.derep} \
      --batch={params.batch} --log={log} --config={params.config}"""
    
rule remove_chimeras:
  input:
    asv = expand("output/dada2/asv_batch/{batch}_asv.qs", batch = batches)
  output:
    asv = "output/dada2/remove_chim/asv_mat_wo_chim.qs",
    summary = "output/dada2/remove_chim/asv_mat_wo_chim_summary.qs",
    matspy = "workflow/report/dada2qc/asv_matrix_wo_chim.png"
  threads: config["threads"]
  params:
    config = "config/config.yaml"
  log:
    "logs/dada2/05_remove_chimeras.log"
  shell:
    """Rscript workflow/scripts/dada2/remove_chimeras.R \
    {output.asv} {output.summary} {output.matspy} \
    {input.asv} \
    --log={log} --config={params.config} --cores={threads}"""

# rule filter_asvs:
#   input:
#     seqtab = rules.remove_chimeras.output.asvs,
#     negcontrol = config["negcontroltable"]
#   output:
#     plot_seqlength = "figures/qc/nasvs_by_seqlength.png",
#     plot_seqabundance = "figures/qc/nasvs_by_seqabundance.png",
#     seqtab_filt = "data/asv/seqtab_nochimeras_qc.qs"
#   log:
#     "logs/dada2/filter_qc.txt"
#   script:
#     "../scripts/dada2/filter_asvs.R"

# rule stats:
#   input:
#     nreads_filtered = "data/stats/Nreads_filtered.txt",
#     nreads_dereplicated = "data/stats/Nreads_dereplicated.txt",
#     nreads_chim_removed = "data/stats/Nreads_chimera_removed.txt"
#   output:
#     nreads = "data/stats/Nreads_dada2.txt",
#     fig_step = "figures/qc/dada2steps_vs_abundance.png",
#     fig_step_rel = "figures/qc/dada2steps_vs_relabundance.png",
#   log:
#     "logs/dada2/summarize_stats.txt"
#   script:
#     "../scripts/dada2/summarize_nreads.R"

# # rule filterLength:
# #     input:
# #         seqtab= rules.removeChimeras.output.rds
# #     output:
# #         plot_seqlength= "figures/Lengths/Sequence_Length_distribution.pdf",
# #         plot_seqabundance= "figures/Lengths/Sequence_Length_distribution_abundance.pdf",
# #         rds= "output/seqtab.rds",
# #         tsv=  "output/seqtab.tsv",
# #     threads:
# #         1
# #     conda:
# #         "../envs/dada2.yaml"
# #     log:
# #         "logs/dada2/filterLength.txt"
# #     script:
# #         "../scripts/dada2/filterLength.R"



# # rule IDtaxa:
# #     input:
# #         seqtab= "output/seqtab.rds",
# #         ref= lambda wc: config['idtaxa_dbs'][wc.ref]
# #     output:
# #         taxonomy= "taxonomy/{ref}.tsv",
# #     threads:
# #         config["threads"]
# #     log:
# #         "logs/dada2/IDtaxa_{ref}.txt"
# #     script:
# #         "../scripts/dada2/IDtaxa.R"

# # localrules: get_ggtaxonomy

# # rule get_ggtaxonomy:
# #     input:
# #         rules.IDtaxa.output
# #     output:
# #         taxonomy= "taxonomy/{ref}_gg.tsv",
# #     threads:
# #         1
# #     log:
# #         "logs/dada2/get_ggtaxonomy_{ref}.txt"
# #     run:

# #         import pandas as pd

# #         tax = pd.read_csv(input[0],sep='\t',index_col=0)

# #         out= tax.dropna(how='all',axis=1).dropna(how='all',axis=0)

# #         out= out.apply(lambda col: col.name[0]+'_'+col.dropna())
# #         out= out.apply(lambda row: ';'.join(row.dropna()),axis=1)
# #         out.name='Taxonomy'


# #         out.to_csv(output[0],sep='\t',header=True)







# # rule get_rep_seq:
# #     input:
# #         "output/seqtab.tsv",
# #     output:
# #         "taxonomy/rep_seq.fasta"
# #     run:
# #         with open(input[0]) as infile:
# #             seqs = infile.readline().strip().replace('"','').split()
# #         with open(output[0],'w') as outfile:
# #             for s in seqs:
# #                 outfile.write(f">{s}\n{s}\n")
