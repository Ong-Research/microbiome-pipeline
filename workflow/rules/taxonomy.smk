rule extract_fasta:
  input:
    asv = "output/dada2/remove_chim/asv_mat_wo_chim.qs"
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
  resources:
    mem_mb = 100000
  shell:
    """kraken2 --db {input.ref} --threads {threads} \
      --output {output.out} --use-names \
      --report {output.summary} --use-mpa-style \
      --confidence {params.confidence} \
      --classified-out {output.classified} \
      --memory-mapping \
      --unclassified-out {output.unclassified} {input.fasta}"""

rule blast_fast_single:
  input:
    fasta = "output/taxa/fasta/asv_sequences.fa"
  output:
    blast = "output/taxa/blast/{db}/blast_results.tsv"
  params:
    db = "data/blast_dbs/{db}",
    fmt = config["blast"]["format"],
    perc = config["blast"]["perc"]
  threads: 24
  shell:
    """blastn -db {params.db} -query {input.fasta} \
      -out {output.blast} -outfmt {params.fmt} \
      -perc_identity {params.perc} \
      -num_threads {threads}"""


rule clean_kraken_with_blast:
  input:
    kraken = "output/taxa/kraken/{ref}/kraken_results.out",
    kraken_summary = "output/taxa/kraken/{ref}/kraken_summary.out",
    blast = "output/taxa/blast/{db}/blast_results.tsv",
    fasta = "output/taxa/fasta/asv_sequences.fa"
  output:
    taxa = "output/taxa/kraken_match/{ref}/{db}/kraken_rdata.qs",
    hits = "output/taxa/kraken_match/{ref}/{db}/kraken_hits.qs"
  threads: 1
  log:
    "logs/taxonomy/{ref}/{db}/clean_kraken2_w_blast.log"
  shell:
    """Rscript workflow/scripts/taxonomy/clean_kraken_wblast.R \
      {output.taxa} {output.hits} \
      --kraken={input.kraken} --summary={input.kraken_summary} \
      --blast={input.blast} --fasta={input.fasta} \
      --log={log} --cores={threads}
    """

rule merge_kraken:
  input:
    taxa = expand("output/taxa/kraken_match/{ref}/{db}/kraken_rdata.qs",
      ref = kraken_dbs, db = blast_dbs),
    hits = expand("output/taxa/kraken_match/{ref}/{db}/kraken_hits.qs",
      ref = kraken_dbs, db = blast_dbs),
    blast = expand("output/taxa/blast/{db}/blast_results.tsv", db = blast_dbs)
  output:
    taxa = "output/taxa/kraken_merged/kraken_rdata.qs",
    hits = "output/taxa/kraken_merged/kraken_hits.qs"
  params:
    config = "config/config.yaml"
  threads: 1
  script:
    """../scripts/taxonomy/merge_dblabels_wblast.R"""


# generates the taxa ID -> name mapping file
# as some databases do not use NCBI
# rule parse_kraken_summary:
#   input:
#     kraken = "output/taxa/kraken/{ref}/kraken_summary.out"
#   output:
#     standard = "output/taxa/kraken/{ref}/kraken_taxmap_standard.tsv",
#     full = "output/taxa/kraken/{ref}/kraken_taxmap_full.tsv"
#   threads: 1
#   log:
#     "logs/taxonomy/parse_kraken2_summary_{ref}.txt"
#   shell:
#     """Rscript workflow/scripts/taxonomy/parse_kraken_summary.R \
#       {output.standard} {output.full} {input.kraken} \
#      --log={log}"""

# rule parse_kraken:
#   input:
#     kraken = "output/taxa/kraken/{ref}/kraken_results.out",
#     mapfile = "output/taxa/kraken/{ref}/kraken_taxmap_standard.tsv"
#   output:
#     taxa = "output/taxa/kraken/{ref}/kraken_taxatable.qs",
#     summary = "output/taxa/kraken/{ref}/kraken_taxasummary.tsv"
#   threads: config["threads"] / 2
#   log:
#     "logs/taxonomy/parse_kraken2_{ref}.txt"
#   shell:
#     """Rscript workflow/scripts/taxonomy/parse_kraken.R \
#       {output.taxa} {output.summary} {input.kraken} \
#       --map {input.mapfile} \
#       --log={log} --cores={threads}"""

