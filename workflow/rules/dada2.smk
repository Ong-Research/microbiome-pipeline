rule filter_and_trim_sample:
  input:
    end1 = lambda wc: sample_dict[wc.sample]["end1"],
    end2 = lambda wc: sample_dict[wc.sample]["end2"]
  output:
    end1 = temp("output/dada2/filtered/{sample}_filtered_R1.fastq.gz"),
    end2 = temp("output/dada2/filtered/{sample}_filtered_R2.fastq.gz"),
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
      --batch={params.batch} --log={log} --config={params.config}
      
      if [[ -s {output.summary} && ! -s {output.end1} && ! -s {output.end2} ]]; then
        echo "Filtering succeeded but all reads were removed. Creating temp placeholder files."
        touch {output.end1}
        touch {output.end2}
      fi
      """

rule learn_error_rates_batch_end:
  input:
    filtered = lambda wc: expand("output/dada2/filtered/{sample}_filtered_{end}.fastq.gz", sample = get_batch_samples(sample_table, wc.batch), end = wc.end)
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
    merge = temp("output/dada2/merge/{batch}/{sample}_asv.qs")
  params:
    config = "config/config.yaml",
    sample_name = lambda wc: wc.sample,
    batch = lambda wc: sample_dict[wc.sample]["batch"]
  threads: 1
  log:
    "logs/dada2/03_{batch}_{sample}_merge.log"
  shell:
    """
      if [[ -e {input.filt_end1} && -e {input.filt_end2} && ! -s {input.filt_end1} && ! -s {input.filt_end2} ]]; then
          touch {output.merge} 
          echo "Filtering succeeded, but no reads were left after filtering:"
          ls -l {input.filt_end1} {input.filt_end2}
          echo "Generating placeholder output file: {output.merge}"
      else 
        Rscript workflow/scripts/dada2/dereplicate_one_sample_pair.R \
        {output.merge} {params.sample_name} \
        {input.filt_end1} {input.filt_end2} \
        --end1_err={input.mat_end1} --end2_err={input.mat_end2} \
        --log={log} --batch={params.batch} --config={params.config}
      fi
      """

rule dereplicate_batch:
  input:
    derep = lambda wc: expand("output/dada2/merge/{batch}/{sample}_asv.qs", sample = get_batch_samples(sample_table, wc.batch), batch = wc.batch)
  output:
    asv = "output/dada2/asv_batch/{batch}_asv.qs",
    summary = "output/dada2/asv_batch/{batch}_summary.tsv"
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

rule filter_asvs:
  input:
    seqtab = "output/dada2/remove_chim/asv_mat_wo_chim.qs",
    negcontrol = config["negcontroltable"] # "data/negcontrols.qs"
  output:
    plot_seqlength = "workflow/report/dada2qc/nasvs_by_seqlength.png",
    plot_seqabundance = "workflow/report/dada2qc/nasvs_by_seqabundance.png",
    seqtab_filt = "output/dada2/after_qc/asv_mat_wo_chim.qs"
  params:
    config = "config/config.yaml"
  log:
    "logs/dada2/06_filter_qc.log"
  threads: 1
  shell:
    """Rscript workflow/scripts/dada2/filter_asvs.R \
      {output.seqtab_filt} \
      {output.plot_seqlength} {output.plot_seqabundance} \
      {input.seqtab} {input.negcontrol} \
      --log={log} --config={params.config} --cores={threads}"""

rule collect_filtered:
  input:
    filt = expand("output/dada2/summary/{sample}_summary_filtered.tsv", sample = all_samples)
  output:
    filt = "output/dada2/filtered/all_sample_summary.tsv"
  log:
    "logs/dada2/07_collect_filtered.log"
  threads: 1
  shell:
    """Rscript workflow/scripts/dada2/collect_tsv_files.R \
      {output.filt} {input.filt} \
      --log={log}"""

rule collect_merged:
  input:
    merg = expand("output/dada2/asv_batch/{batch}_summary.tsv",
      batch = batches)
  output:
    merg = "output/dada2/asv_batch/all_sample_summary.tsv"
  log:
    "logs/dada2/07_collect_derep.log"
  threads: 1
  shell:
    """Rscript workflow/scripts/dada2/collect_tsv_files.R \
      {output.merg} {input.merg} \
      --log={log}"""

rule collect_summary:
  input:
    filt = "output/dada2/filtered/all_sample_summary.tsv",
    merg = "output/dada2/asv_batch/all_sample_summary.tsv",
    asvs = "output/dada2/after_qc/asv_mat_wo_chim.qs"
  output:
    nreads = "output/dada2/stats/Nreads_dada2.txt",
    fig_step = "workflow/report/dada2qc/dada2steps_vs_abundance.png",
    fig_step_rel = "workflow/report/dada2qc/dada2steps_vs_relabundance.png"
  log:
    "logs/dada2/08_summarize_stats.txt"
  shell:
    """Rscript workflow/scripts/dada2/summarize_nreads.R \
      {output.nreads} {output.fig_step} {output.fig_step_rel} \
      {input.filt} {input.merg} {input.asvs} \
      --log={log}"""
