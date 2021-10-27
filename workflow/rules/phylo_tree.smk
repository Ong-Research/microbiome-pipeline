rule import_fasta:
  input:
    fa = "output/taxa/fasta/asv_sequences.fa"
  output:
    qza = "output/phylotree/asv_sequences.qza"
  conda:
    "../envs/qiime2_2021_08_conda.yaml"
  shell:
    """qiime tools import \
      --input-path {input.fa} \
      --output-path {output.qza} \
      --type 'FeatureData[Sequence]'"""
    
rule make_tree:
  input:
    qza = "output/phylotree/asv_sequences.qza"
  output:
    alignment = "output/phylotree/aligned_sequences.qza",
    mask_alignment = "output/phylotree/masked_aligned_sequences.qza",
    unroot_tree = "output/phylotree/unrooted_tree.qza",
    root_tree = "output/phylotree/rooted_tree.qza"
  conda:
    "../envs/qiime2_2021_08_conda.yaml"
  shell:
    """qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences {input.qza} \
  --o-alignment {output.alignment} \
  --o-masked-alignment {output.mask_alignment} \
  --o-tree {output.unroot_tree} \
  --o-rooted-tree {output.root_tree}"""

rule export_tree:
  input:
    root_tree = "output/phylotree/rooted_tree.qza"
  output:
    newick_tree = "output/phylotree/newick/tree.nwk"
  conda:
    "../envs/qiime2_2021_08_conda.yaml"
  shell:
    """qiime tools export \
      --input-path {input.root_tree} \
      --output-path output/phylotree/newick
    """
    # this will make a directory newick/tree.nwk that contains
    # a file tree.nwk... we don't need that
      #--output-path {output.newick_tree}"""
