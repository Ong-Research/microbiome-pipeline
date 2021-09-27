#!/usr/local/bin/Rscript

#' Prepara mia object
#' @author rwelch2

"Prepare mia object

Usage:
prepare_mia_object.R [<mia_file>] [--asv=<asv_file> --taxa=<taxa_file> --tree=<tree_file> --meta=<meta_file>] [--asv_prefix=<prefix> --log=<logfile>]
prepare_mia_object.R (-h|--help)
prepare_mia_object.R --version

Options:
-h --help    show this screen
--asv=<asv_file>    ASV matrix file
--taxa=<taxa_file>    Taxa file
--tree=<tree_file>    Tree file
--meta=<meta_file>    Metadata file
--log=<logfile>    name of the log file [default: logs/filter_and_trim.log]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "prepare mia file V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$asv <- "output/dada2/after_qc/asv_mat_wo_chim.qs"
  arguments$taxa <- "output/taxa/kraken/minikraken/kraken_taxatable.qs"
  arguments$tree <- "output/phylotree/newick/tree.nwk"
  arguments$meta <- "data/meta.tsv"
  arguments$asv_prefix <- "HSD2M"
}

print(arguments)

info <- Sys.info();

message("loading packages")
library(magrittr)
library(tidyverse)
library(TreeSummarizedExperiment)
library(Biostrings)
library(ape)
library(mia)
library(qs)

stopifnot(
  file.exists(arguments$asv),
  file.exists(arguments$taxa),
  file.exists(arguments$tree),
  file.exists(arguments$meta)
)

asv <- qs::qread(arguments$asv)
taxa <- qs::qread(arguments$taxa)
tree <- ape::read.tree(arguments$tree)
meta <- readr::read_tsv(arguments$meta)

asv_sequences <- colnames(asv)
colnames(asv) <- str_c(arguments$asv_prefix, seq_along(asv_sequences),
  sep = "_")
names(asv_sequences) <- colnames(asv)
asv_sequences <- Biostrings::DNAStringSet(asv_sequences)

asv_aux <- tibble::tibble(asv = colnames(asv))


# need to fix parse_taxa to return a tibble of length 
# equal to all the sequences and not only the sequences with known information

cdata <- meta %>%
  as.data.frame() %>%
  tibble::column_to_rownames("key")

out <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(Count = t(asv)),
  colData = cdata,
  rowData = asv_aux %>%
    dplyr::left_join(taxa, by = "asv") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("asv"),
  rowTree = tree)

metadata(out)[["date_processed"]] <- Sys.Date()
metadata(out)[["sequences"]] <- asv_sequences

fs::dir_create(dirname(arguments$mia_file))
qs::qsave(out, arguments$mia_file)
