#!/usr/local/bin/Rscript

#' Parses kraken2 results into a readable table.
#' By default, maps tax IDs to taxonomy strings using the NCBI database.
#' If --map file provided (output of parse_kraken_summary), uses that as a mapping instead.
#' Output files are taxa_table and taxa_summary
#' @param kraken results computed by kraken2
#' @author rwelch
#' @author chasman

"Parse kraken2 results

Usage:
parse_kraken.R [<taxa_table> <taxa_summary>] [<kraken_file>] [--map=<map_file> --log=<logfile> --cores=<cores>]
parse_kraken.R (-h|--help)
parse_kraken.R --version

Options:
--log=<logfile>    name of the log file [default: ./parse_kraken.log]
--map=<map_file>   name of file with map from taxid to taxonomy [default: NULL]
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

  db = "silva"
  arguments$taxa_table <- "parse_output"
  arguments$taxa_summary <- "parse_output_summary"
  arguments$kraken_file <- sprintf("output/taxa/kraken/%s/kraken_results.out", db)
  arguments$map <- "output_test"
}


info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(taxizedb)
library(BiocParallel)

bpp <- BiocParallel::MulticoreParam(workers = as.numeric(arguments$cores))

tax_order <- c("domain","superkingdom", "kingdom", "phylum", 
    "class", "order", "family", "genus", "species")

# mapping file provided
if (!is.null(arguments$map)) {

  message("reading labels from mapping file ", arguments$map)

  # taxonomy info
  map_tb <- readr::read_tsv(arguments$map)
  # ASV assignments
  kraken <- readr::read_tsv(arguments$kraken_file, 
    col_names = c("status", "asv", "taxid", "seq_len", "lca_mapping"))
  labels <- kraken %>%
    dplyr::select(asv, taxid) %>%
    dplyr::left_join(map_tb) %>%
    dplyr::select(-taxname, -taxid)

} else {
  message("getting ncbi database")
  ncbi_db <- taxizedb::db_download_ncbi()

  message("reading labels from kraken")
  kraken <- readr::read_tsv(arguments$kraken_file, 
    col_names = c("status", "asv", "id", "seq_length", "id_bp"))

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

    # any_of will allow for missing taxonomic levels (eg domain)
    labels %>%
      dplyr::select(asv, id_taxa_clean) %>%
      tidyr::unnest(cols = c(id_taxa_clean)) %>%
      dplyr::select(asv, tidyselect::any_of(taxa)) 

  }

  message("parsing labels")
  labels <- get_taxa_from_id(kraken) 
  labels <- clean_id_taxa_wrap(labels, tax_order, bpp)  
  # after the cleaning, we have a table like this:
  #    asv    phylum         class           order             family  genus species
  #   <chr>  <chr>          <chr>           <chr>             <chr>   <chr> <chr>  
  #  1 asv_1  Chordata       Mammalia        Artiodactyla      Cervid… Cerv… Cervus…
  #  2 asv_2  Chordata       Mammalia        Artiodactyla      Cervid… Cerv… Cervus…
  #  3 asv_3  Actinobacteria Actinomycetia   Pseudonocardiales Pseudo… Sacc… Saccha…
}

taxa_summary <- labels %>%
  tidyr::pivot_longer(tidyselect::any_of(tax_order), 
    names_to = "taxa", values_to = "value") %>%
  dplyr::mutate(
    taxa = factor(taxa, levels = tax_order)) %>%
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

