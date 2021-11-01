#!/usr/local/bin/Rscript

#' Parses kraken2 results into a readable table
#' @param kraken results computed by kraken2
#' @author rwelch

"Parse kraken2 results

Usage:
parase_kraken.R [<taxa_table> <taxa_summary>] [<kraken_file>] [--log=<logfile> --cores=<cores>]
parase_kraken.R (-h|--help)
parase_kraken.R --version

Options:
--log=<logfile>    name of the log file [default: ./parse_kraken.log]
--cores=<cores>    number of parallel CPUs [default: 8]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "parse kraken2 results V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$taxa_table <- "output"
  #arguments$kraken_file <- "output/taxa/kraken/minikraken/kraken_results.out"
  arguments$kraken_file <- "output/taxa/kraken/silva/kraken_results.out"

}


info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(taxizedb)
library(BiocParallel)

bpp <- BiocParallel::MulticoreParam(workers = as.numeric(arguments$cores))

message("getting ncbi database")
ncbi_db <- taxizedb::db_download_ncbi()
taxa <- c("phylum", "class", "order", "family", "genus", "species")

message("reading labels from kraken")
kraken <- readr::read_tsv(arguments$kraken_file, col_names = FALSE)

get_taxa_from_id <- function(results) {
  query <- taxizedb::classification(results$id, db = "ncbi")

  results %>%
    dplyr::mutate(
      id_taxa = map(query, list),
      id_taxa = purrr::map(id_taxa, ~ unique(.[[1]])),
      is_na = ! purrr::map_lgl(id_taxa, is.data.frame))

}

clean_id_taxa <- function(id_taxa, taxa) {

  name <- NULL

  if (is.data.frame(id_taxa)) {
    id_taxa %<>%
      dplyr::filter(rank %in% taxa) %>%
      dplyr::select(-id)
    id_taxa %<>%
      tidyr::pivot_wider(names_from = rank, values_from = name)
  }
  id_taxa
}

clean_id_taxa_wrap <- function(labels, taxa, bpp) {

  labels %<>%
    dplyr::mutate(
      id_taxa_clean = BiocParallel::bplapply(id_taxa, clean_id_taxa, taxa,
        BPPARAM = bpp))

  labels %>%
    dplyr::select(asv, id_taxa_clean) %>%
    tidyr::unnest(cols = c(id_taxa_clean)) %>%
    dplyr::select(asv, tidyselect::one_of(taxa))

}

message("parsing labels")
kraken %<>%
  rlang::set_names(c("rank", "asv", "id", "seq_length", "id_bp"))
  
labels <- get_taxa_from_id(kraken)
labels <- clean_id_taxa_wrap(labels, taxa, bpp)

taxa_summary <- labels %>%
  tidyr::pivot_longer(-asv, names_to = "taxa", values_to = "value") %>%
  dplyr::mutate(
    taxa = factor(taxa, levels = c("phylum", "class", "order", "family",
      "genus", "species"))) %>%
  dplyr::group_by(taxa) %>%
  dplyr::summarize(
    total = length(value),
    label = sum(!is.na(value)), .groups = "drop") %>%
  dplyr::mutate(perc = label / total)

message("saving results")

fs::dir_create(dirname(arguments$taxa_table))
fs::dir_create(dirname(arguments$taxa_summary))

labels %>%
  qs::qsave(arguments$taxa_table)
  
taxa_summary %>%
  readr::write_tsv(arguments$taxa_summary)
