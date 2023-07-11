
rule plot_quality_profiles:
  input:
    end1 = lambda wc: sample_dict[wc.sample]["end1"],
    end2 = lambda wc: sample_dict[wc.sample]["end2"]
  log:
    "logs/qc/plot_qc_profiles_{sample}.log"
  output:
    "workflow/report/quality_profiles/{sample}.png"
  shell:
    """Rscript workflow/scripts/dada2/plot_quality_profiles.R \
      {output} --end1={input.end1} --end2={input.end2} \
      --log={log}"""

rule fastqc:
  input:
    R1 = sample_table.end1.values,
    R2 = sample_table.end2.values
  params:
    threads = config["threads"]
  conda:
    "../envs/fastqc.yaml"
  output:
    zip=temp(expand("output/quality_control/fastqc/{sample}_fastqc.zip", sample = all_files))
  log:
    "logs/qc/fastqc.txt"
  shell:
    """fastqc -o output/quality_control/fastqc -t {params.threads} {input.R1} {input.R2}"""
    
rule multiqc:
  input:
    expand("output/quality_control/fastqc/{sample}_fastqc.zip", sample = all_files)
  output:
    "output/quality_control/multiqc/multiqc_report.html"
  conda:
    "../envs/multiqc.yaml"
  shell:
    """multiqc output/quality_control/fastqc -o output/quality_control/multiqc"""