#!/usr/local/bin/Rscript

#' Merges a list of kraken2 labels based on different dbs using auxiliary blast
#' alignments and hits
#'
#' Output files are taxa_table and taxa_summary
#' @param kraken results computed by kraken2
#' @author rwelch


# save(snakemake, file = "zz.RData")
# load("zz.RData")


stopifnot(
  all(file.exists(snakemake@input[["taxa"]])),
  all(file.exists(snakemake@input[["hits"]])),
  all(file.exists(snakemake@input[["blast"]])))

library(magrittr)
library(tidyverse)
library(vroom)
library(microbiome.misc)
library(progress)

# snakemake
rdata_files <- snakemake@input[["taxa"]]
hits_files <- snakemake@input[["hits"]]
blast_files <- snakemake@input[["blast"]]

config <- here::here(snakemake@params[["config"]]) %>%
  yaml::read_yaml()


blast_dbs <- basename(dirname(blast_files))

format <- stringr::str_split(config[["blast"]][["format"]], " ")
format <- format[[1]][-1]
format <- stringr::str_remove_all(format, "\\'")
format[1] <- "asv"

blast_list <- tibble::tibble(blast_db = blast_dbs,
  blast = purrr::map(blast_files, vroom::vroom,
  col_names = format))

kraken_blast_min_perc <- config[["kraken_blast_min_perc"]] * 100

blast_list %<>%
  dplyr::mutate(
    blast = purrr::map(blast, dplyr::filter, pident >= kraken_blast_min_perc),
    blast = purrr::map(blast, tidyr::nest, blast = -c(asv)))

get_db <- microbiome.misc:::get_db

all_hits <- tibble::tibble(
  blast_db = get_db(hits_files, "blast"),
  kraken_db = get_db(hits_files, "kraken"),
  hits = purrr::map(hits_files, qs::qread))

all_labels <- tibble::tibble(
  blast_db = get_db(rdata_files, "blast"),
  kraken_db = get_db(rdata_files, "kraken"),
  hits = purrr::map(rdata_files, qs::qread))

blast_list %<>%
  tidyr::unnest(cols = c(blast))

db_hierarchy <- config[["kraken_db_merge"]]

message("subsetting kraken dbs: ",
  stringr::str_c(db_hierarchy, collapse = ", "))

all_hits %<>%
  dplyr::filter(kraken_db %in% db_hierarchy) %>%
  tidyr::unnest(cols = c(hits)) %>%
  tidyr::nest(hits = -c(blast_db, asv))

all_labels %<>%
  dplyr::filter(kraken_db %in% db_hierarchy) %>%
  tidyr::unnest(cols = c(hits)) %>%
  tidyr::nest(labels = -c(blast_db, asv))

merge_labels <- purrr::reduce(
  list(all_hits, all_labels, blast_list), dplyr::full_join,
    by = c("blast_db", "asv"))

pb <- progress::progress_bar$new(total = nrow(merge_labels))

match_asvs_blast_pbar <- function(hits, labels, blast, hierarchy) {
  pb$tick()
  microbiome.misc::match_asvs_blast(hits, labels, blast, hierarchy)

}


merge_labels %<>%
  dplyr::mutate(
    merged = purrr::pmap(list(hits, labels, blast),
      match_asvs_blast_pbar, config[["kraken_db_merge"]]))

if (length(unique(merge_labels$blast_db)) > 1) {
  stop("more than 1 db, fix for next time")
}

merge_labels %<>%
  dplyr::mutate(
    hits_merged = purrr::map(merged, "hits"),
    labels_merged = purrr::map(merged, "labels"))

merge_hits <- merge_labels %>%
  dplyr::select(asv, hits_merged) %>%
  tidyr::unnest(cols = c(hits_merged)) %>%
  dplyr::select(kraken_db, tidyselect::everything())

merge_hits %>%
  qs::qsave(snakemake@output[["hits"]])

merge_lbs <- merge_labels %>%
  dplyr::select(asv, labels_merged) %>%
  tidyr::unnest(cols = c(labels_merged)) %>%
  dplyr::select(kraken_db, tidyselect::everything())

merge_lbs %>%
  qs::qsave(snakemake@output[["taxa"]])
