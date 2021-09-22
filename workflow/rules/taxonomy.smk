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

# rule kraken_taxonomy:
#   input:
#     fasta = rules.extract_fasta.output.fasta,
#     ref = lambda wc: config["kraken_dbs"][wc.ref]
#   output:
#     out = "data/taxonomy/kraken_{ref}.out",
#     summary = "data/taxonomy/kraken_{ref}.summary",
#     classified = temp("data/taxonomy/kraken_{ref}.classified"),
#     unclassified = temp("data/taxonomy/kraken_{ref}.unclassified")
#   threads:
#     config["threads"]
#   params:
#     confidence = config["confidence"]
#   log:
#     "logs/taxonomy/kraken2_{ref}.txt"
#   shell:
#     """
#     kraken2 --db {input.ref} --threads {threads} \
#       --output {output.out} --report {output.summary} \
#       --confidence {params.confidence} \
#       --classified-out {output.classified} \
#       --unclassified-out {output.unclassified} {input.fasta}
#     """

# localrules:
#   parse_kraken

# rule parse_kraken:
#   input:
#     kraken = rules.kraken_taxonomy.output.out
#   output:
#     taxa_qs = "data/taxonomy/kraken_{ref}_labels.qs",
#     taxa_tsv = "data/taxonomy/kraken_{ref}_labels.tsv",
#     summary = "data/taxonomy/kraken_{ref}_summary.tsv"
#   threads:
#     config["threads"]
#   log:
#     "logs/taxonomy/parse_kraken2_{ref}.txt"
#   script:
#     "../scripts/taxonomy/parse_kraken.R"
