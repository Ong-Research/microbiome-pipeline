
rule create_tree:
  """
  Computes a phylogenetic tree using the fastreeR package
  """
  input:
    fasta = "output/taxa/fasta/asv_sequences.fa"
  output:
    tree = "output/phylotree/newick/tree.nwk"
  script:
    """../scripts/phylo/compute_tree.R"""
