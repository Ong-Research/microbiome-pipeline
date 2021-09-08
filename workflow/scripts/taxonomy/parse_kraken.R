
source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(taxizedb)
library(furrr)
library(future)

future::plan(multiprocess, workers = snakemake@threads)

message("getting ncbi database")
ncbi_db <- taxizedb::db_download_ncbi()
taxa <- c("phylum", "class", "order", "family", "genus", "species")

message("reading labels from kraken")
kraken <- snakemake@input[["kraken"]]
kraken <- readr::read_tsv(kraken, col_names = FALSE)

get_taxa_from_id <- function(results) {
  query <- taxizedb::classification(results$id, db = "ncbi")

  results %>%
    dplyr::mutate(
      id_taxa = map(query, list),
      id_taxa = purrr::map(id_taxa, ~ unique(.[[1]])),
      is_na = ! purrr::map_lgl(id_taxa, is.data.frame))

}

clean_id_taxa <- function(id_taxa, taxa) {

  `%<>%` <- magrittr::`%<>%`
  name <- NULL

  if (is.data.frame(id_taxa)) {
    id_taxa %<>%
      dplyr::filter(rank %in% taxa) %>%
      dplyr::select(-id)
    id_taxa %<>% tidyr::spread(rank, name)
	}
  id_taxa
}

clean_id_taxa_wrap <- function(labels, taxa) {

  labels %<>%
    dplyr::mutate(
      id_taxa_clean = furrr::future_map(id_taxa, clean_id_taxa, taxa))

  labels %>%
    dplyr::select(asv, id_taxa_clean) %>%
    tidyr::unnest(cols = c(id_taxa_clean)) %>%
    dplyr::select(asv, tidyselect::one_of(taxa))

}


message("parsing labels")
kraken %<>%
  rlang::set_names(c("rank", "asv", "id", "seq_length", "id_bp"))
labels <- get_taxa_from_id(kraken)
labels <- clean_id_taxa_wrap(labels, taxa)

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

labels %>% qs::qsave(snakemake@output[["taxa_qs"]])
labels %>% readr::write_tsv(snakemake@output[["taxa_tsv"]])
taxa_summary %>% readr::write_tsv(snakemake@output[["summary"]])