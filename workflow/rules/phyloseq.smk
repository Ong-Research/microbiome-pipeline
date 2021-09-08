
rule phylo_tree:
  input:
    asv="data/asv/seqtab_nochimeras_qc.qs"
  output:
    alignment="data/fasta/sequence_alignment.fa",
    tree="data/phyloseq/tree.nwk"
  log:
    "logs/phyloseq/compute_phylotree.txt"
  script:
    "../scripts/phyloseq/align_sequences.R"


rule init_phyloseq:
  input:
    asv="data/asv/seqtab_nochimeras_qc.qs",
    taxa="data/taxonomy/kraken_minikraken2_labels.qs",
    meta=config["metadata"]
    # tree="data/phyloseq/tree.nwk"
  output:
    sequences="data/phyloseq/asv_sequences.qs",
    phyloseq="data/phyloseq/asv_phyloseq.qs"
  log:
    "logs/phyloseq/init_phyloseq.txt"
  script:
    "../scripts/phyloseq/init_phyloseq.R"

rule normalize:
  input:
    phyloseq = "data/phyloseq/asv_phyloseq.qs"
  output:
    phyloseq_norm = "data/phyloseq/asv_phyloseq_norm.qs"
  script:
    "../scripts/deseq2/normalize_samples.R"

rule alpha_div:
  input:
    phyloseq="data/phyloseq/asv_phyloseq.qs"
  output:
    alpha="data/phyloseq/div/alphadiv.qs"
  log:
    "logs/phyloseq/alphadiv.txt"
  script:
    "../scripts/phyloseq/compute_alphadiv.R"

rule beta_div:
  input:
    phyloseq="data/phyloseq/asv_phyloseq_norm.qs"
  output:
    distance="data/phyloseq/div/beta_distance.qs",
    pcoa = "data/phyloseq/div/betadiv.qs"
  threads:
    config["threads"]
  log:
    "logs/phyloseq/betadiv.txt"
  script:
    "../scripts/phyloseq/compute_betadiv.R"
