#!/usr/local/bin/Rscript

#' Parses kraken2 results into a readable table.
#' By default, maps tax IDs to taxonomy strings using the NCBI database.
#' If --mpa file provided (output of KrakenTools kreport2mpa.py), uses that as a mapping instead.
#' @param kraken results computed by kraken2
#' @author rwelch

"Parse kraken2 results

Usage:
parse_kraken.R [<taxa_table> <taxa_summary>] [<kraken_file>] [--mpa=<mpa_file> --log=<logfile> --cores=<cores>]
parse_kraken.R (-h|--help)
parse_kraken.R --version

Options:
--log=<logfile>    name of the log file [default: ./parse_kraken.log]
--mpa=<mpa_file>   name of Metaphlan-style taxonomy string to ID map file [default: NULL]
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
  arguments$taxa_table <- "output_test"
  arguments$taxa_summary <- sprintf("output/taxa/kraken/%s/kraken_summary.out", db)
  arguments$kraken_file <- sprintf("output/taxa/kraken/%s/kraken_results.out", db)
  arguments$mpa <- "tmp" #sprintf("output/taxa/kraken/%s/kraken_mpa.tsv", db)
}


info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(taxizedb)
library(BiocParallel)

bpp <- BiocParallel::MulticoreParam(workers = as.numeric(arguments$cores))

if (!is.null(arguments$mpa) && 
  !is.null(arguments$taxa_summary) && 
  !is.null(arguments$taxa_table)) {

  # these are the levels we are interested in
  # there may be others in the output (like G2, S2), 
  # but we are only going this far
  taxa <- list(D = "domain", K="kingdom", P = "phylum", 
    C = "class", O = "order", F = "family", G = "genus", 
    S = "species")

  #' returns a tibble in format
  # A tibble: 1 × 7
  #   kingdom  phylum           class          order  family  genus  tax_string     
  #   <chr>    <chr>            <chr>          <chr>  <chr>   <chr>  <chr>          
  # 1 Bacteria Actinobacteriota Actinobacteria Bifid… Bifido… Bifid… k__Bacteria|p_…
  clean_mpa <- function(tax_string, tax_tib) {
    tax_split <- stringr::str_split(tax_string, pattern = "\\|") %>% .[[1]]
    tibble::tibble(split = tax_split) %>%
      tidyr::separate(split, into = c("rank", "tax_name"), sep = "__", remove = F) %>%
      left_join(tax_tib, by = "rank") %>%
      dplyr::select(tax_name, taxlevel) %>%
      pivot_wider(names_from = taxlevel, values_from = tax_name) %>%
      dplyr::mutate(tax_string = tax_string) 
  }


  message("reading labels from kraken and krakentools output")
  # taxonomy info
  summary <- readr::read_tsv(arguments$taxa_summary, 
    col_names = c("frac_asvs", "n_asvs_below",
       "n_asvs_assigned", "rank", "taxid", "tax_name"))
  # MPA strings for taxonomy
  mpa_tb <- readr::read_tsv(arguments$mpa, 
    col_names = c("tax_string", "n_asvs_below", "taxid"))
  # ASV assignments
  kraken <- readr::read_tsv(arguments$kraken_file, 
    col_names = c("status", "asv", "taxid", "seq_len", "lca_mapping"))
  
  # the mpa file may only have
  # domain, phylum, class, order, family, genus, species.
  #mpa_tb <- 
  mpa_tb <- summary %>%
    #dplyr::filter(rank %in% names(taxa)) %>%
    left_join(mpa_tb, by = c("n_asvs_below", "taxid"))
    
  # possible tax levels, mapped to longer strings when possible
  possible <- purrr::map(unique(mpa_tb$rank), 
    ~ ifelse(.x %in% names(taxa), taxa[.x], .x)) %>% 
    unlist
  labels <- kraken %>%
    left_join(mpa_tb, by = "taxid") %>%
    dplyr::select(asv, rank, tax_name, tax_string) %>%
    tidyr::separate(tax_string, 
      into = possible, 
      extra = "drop",
      sep = "\\|", remove = F) 
      #FFFFFF
  

  # select a subset of ranks
  labels %<>% 
    dplyr::select(asv, rank, tax_name, tax_string, all_of(unlist(taxa)))
  # remove the k__/etc prefixes
  labels %<>%
    dplyr::mutate(across(possible, ~ substr(.x, 4, nchar(.x)))) %>%
    dplyr::select(-tax_name, -rank)

} else {
  taxa <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  message("getting ncbi database")
  ncbi_db <- taxizedb::db_download_ncbi()

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

    # any_of will allow for missing taxonomic levels (eg domain)
    labels %>%
      dplyr::select(asv, id_taxa_clean) %>%
      tidyr::unnest(cols = c(id_taxa_clean)) %>%
      dplyr::select(asv, tidyselect::any_of(taxa)) 

  }

  message("parsing labels")
  kraken %<>%
    rlang::set_names(c("rank", "asv", "id", "seq_length", "id_bp"))
    
  labels <- get_taxa_from_id(kraken)  # uses NCBI
  labels <- clean_id_taxa_wrap(labels, taxa, bpp)  
  # after the cleaning, we have a table like this:
  #    asv    phylum         class           order             family  genus species
  #   <chr>  <chr>          <chr>           <chr>             <chr>   <chr> <chr>  
  #  1 asv_1  Chordata       Mammalia        Artiodactyla      Cervid… Cerv… Cervus…
  #  2 asv_2  Chordata       Mammalia        Artiodactyla      Cervid… Cerv… Cervus…
  #  3 asv_3  Actinobacteria Actinomycetia   Pseudonocardiales Pseudo… Sacc… Saccha…
}

tax_order = taxa
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
