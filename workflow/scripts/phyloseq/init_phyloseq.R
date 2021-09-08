source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

library(magrittr)
library(phyloseq)
library(ape)
library(qs)

### define directories
message("loading data")
asv_table <- qs::qread(snakemake@input[["asv"]])
taxa <- qs::qread(snakemake@input[["taxa"]])
# tree <- ape::read.tree(snakemake@input[["tree"]])

tree <- NULL

message("using asv_table in ", snakemake@input[["asv"]])
message("using taxa in ", snakemake@input[["taxa"]])
message("using tree in ", snakemake@input[["tree"]])

seqs <- tibble::tibble(
  asv = stringr::str_c("asv", seq_len(ncol(asv_table)), sep = "_"),
  seqs = colnames(asv_table))

colnames(asv_table) <- dplyr::pull(seqs, asv)

message("creating phyloseq")

if (is.null(taxa) & is.null(tree)) {

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(asv_table, taxa_are_rows = FALSE))

} else if (is.null(tree)) {

  # keep only labelled
  asv_table <- asv_table[, colnames(asv_table) %in% taxa$asv]
  taxa_mat <- as.data.frame(taxa)
  taxa_mat <- tibble::column_to_rownames(taxa_mat, "asv")
  taxa_mat <- as.matrix(taxa_mat)

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(asv_table, taxa_are_rows = FALSE),
    tax_table(taxa_mat))

} else {

  # keep only labelled
  asv_table <- asv_table[, colnames(asv_table) %in% taxa$asv]
  taxa_mat <- as.data.frame(taxa)
  taxa_mat <- tibble::column_to_rownames(taxa_mat, "asv")
  taxa_mat <- as.matrix(taxa_mat)

  tree <- ape::keep.tip(tree, rownames(taxa_mat))

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(asv_table, taxa_are_rows = FALSE),
    phyloseq::tax_table(taxa_mat),
    phyloseq::phy_tree(tree))

}

meta <- readr::read_tsv(snakemake@input[["meta"]])

meta %<>%
  dplyr::distinct() %>%
  as.data.frame() %>%
  tibble::column_to_rownames("sample")

sample_data(ps) <- meta


message("saving sequences in ", snakemake@output[["sequences"]])
seqs %>% qs::qsave(snakemake@output[["sequences"]])

message("saving phyloseq in ", snakemake@output[["phyloseq"]])
ps %>% qs::qsave(snakemake@output[["phyloseq"]])
