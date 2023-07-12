
library(fastreeR)
library(ape)

fasta <- snakemake@input[["fasta"]]

message("computes the distance matrix between fasta sequences")
fasta_dist <- fastreeR::fasta2dist(fasta, kmer = 6)

message("converting distance to tree")
tree <- dist2tree(fasta_dist)
tree <- ape::read.tree(text = tree)

write.tree(tree, file = snakemake@output[["tree"]])
