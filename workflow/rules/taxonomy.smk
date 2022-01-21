rule extract_fasta:
  input:
    asv = "output/dada2/after_qc/asv_mat_wo_chim.qs"
  output:
    fasta = "output/taxa/fasta/asv_sequences.fa"
  params:
    prefix = config["asv_prefix"]
  log:
    "logs/taxonomy/extract_fasta.txt"
  shell:
    """Rscript workflow/scripts/taxonomy/extract_fasta.R \
      {output.fasta} {input.asv} \
      --prefix={params.prefix} --log={log}"""

rule kraken_taxonomy:
  input:
    fasta = "output/taxa/fasta/asv_sequences.fa",
    ref = lambda wc: config["kraken_dbs"][wc.ref]
  output:
    out = "output/taxa/kraken/{ref}/kraken_results.out",
    summary = "output/taxa/kraken/{ref}/kraken_summary.out",
    classified = "output/taxa/kraken/{ref}/kraken_asvs.classified",
    unclassified = "output/taxa/kraken/{ref}/kraken_asvs.unclassified"
  threads:
    config["threads"]
  conda:
    "../envs/kraken.yaml"
  params:
    confidence = config["kraken_confidence"]
  log:
    "logs/taxonomy/kraken2_{ref}.txt"
  shell:
    """kraken2 --db {input.ref} --threads {threads} \
      --output {output.out} --report {output.summary} \
      --confidence {params.confidence} \
      --classified-out {output.classified} \
      --unclassified-out {output.unclassified} {input.fasta}"""

# generates the taxa ID -> name mapping file
# as some databases do not use NCBI
rule parse_kraken_summary:
  input:
    kraken = "output/taxa/kraken/{ref}/kraken_summary.out"
  output:
    standard = "output/taxa/kraken/{ref}/kraken_taxmap_standard.tsv",
    full = "output/taxa/kraken/{ref}/kraken_taxmap_full.tsv"
  threads: 1
  log:
    "logs/taxonomy/parse_kraken2_summary_{ref}.txt"
  shell:
    """Rscript workflow/scripts/taxonomy/parse_kraken_summary.R \
      {output.standard} {output.full} {input.kraken} \
     --log={log}"""

rule parse_kraken:
  input:
    kraken = "output/taxa/kraken/{ref}/kraken_results.out",
    mapfile = "output/taxa/kraken/{ref}/kraken_taxmap_standard.tsv"
  output:
    taxa = "output/taxa/kraken/{ref}/kraken_taxatable.qs",
    summary = "output/taxa/kraken/{ref}/kraken_taxasummary.tsv"
  threads: config["threads"] / 2
  log:
    "logs/taxonomy/parse_kraken2_{ref}.txt"
  shell:
    """Rscript workflow/scripts/taxonomy/parse_kraken.R \
      {output.taxa} {output.summary} {input.kraken} \
      --map {input.mapfile} \
      --log={log} --cores={threads}"""

