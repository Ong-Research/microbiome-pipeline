#!/usr/local/bin/Rscript

#' Cleans kraken2 result into a clean table to be used for `rowData`
#' 
#' Output files are taxa_table and taxa_summary
#' @param kraken results computed by kraken2
#' @author rwelch
#' @author chasman

"Clean kraken2 results with blast

Usage:
clean_kraken_wblast.R [<taxa_table> <hits_table>] [--kraken=<kraken_file> --summary=<summary_file> --blast=<blast_file> --fasta=<fasta_file>] [--log=<logfile> --cores=<cores>]
clean_kraken_wblast.R (-h|--help)
clean_kraken_wblast.R --version

Options:
--log=<logfile>    name of the log file [default: ./parse_kraken.log]
--cores=<cores>    number of parallel CPUs [default: 8]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "clean kraken2 results using blast V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$kraken <- here::here("output/taxa/kraken/k2_std/kraken_results.out")
  arguments$summary <- here::here(
      "output/taxa/kraken/k2_std/kraken_summary.out")
  arguments$blast <- here::here("output/taxa/blast/16S_ribosomal_RNA/blast_results.tsv")
  arguments$fasta <- here::here("output/taxa/fasta//asv_sequences.fa")

}



info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
# library(BiocParallel)

# bpp <- BiocParallel::MulticoreParam(workers = as.numeric(arguments$cores))

tax_order <- c("domain", "phylum",
    "class", "order", "family", "genus", "species")

kraken <- vroom::vroom(arguments$kraken,
  col_names = c("classified", "asv", "taxa", "length", "map"))

ksummary <- vroom::vroom(arguments$summary, col_names = c("taxa", "taxid")) %>%
  dplyr::select(-taxid) %>%
  tidyr::separate(taxa, into = tax_order, sep = "\\|", remove = FALSE,
    fill = "right")

blast <- vroom::vroom(arguments$blast,
  col_names = c("qseqid", "sseqid", "evalue", "bitscore",
    "score", "mismatch", "positive", "stitle", "qframe",
    "sframe", "length", "pident"))

# first remove sequences with width outside fasta sequences
fa <- ShortRead::readFasta(arguments$fasta)
width_vec <- BiocGenerics::width(fa)
width_range <- range(width_vec)

blast %<>%
  dplyr::filter(length >= width_range[1] & length <= width_range[2])

med_mismatches <- median(blast$mismatch)

blast %<>%
  dplyr::filter(mismatch <= med_mismatches) %>%
  dplyr::select(-qframe, -sframe)

kraken_taxids <- kraken %>%
  dplyr::count(taxa) %>%
  dplyr::mutate(
    taxid = stringr::str_split(taxa, "taxid"),
    only_taxa = purrr::map_chr(taxid, 1),
    only_taxa = stringr::str_remove(only_taxa, regex("\\(")),
    only_taxa = stringr::str_trim(only_taxa),
    taxid = purrr::map_chr(taxid, 2),
    taxid = stringr::str_trim(taxid),
    taxid = stringr::str_remove(taxid, "\\)"),
    taxid = as.numeric(taxid))

kraken %<>%
  dplyr::inner_join(
    dplyr::select(kraken_taxids, taxa, taxid), by = "taxa")

# map_to_tibble <- function(map) {

#   map %<>% stringr::str_split("\\:")
#   tibble::tibble(
#     taxid = as.numeric(purrr::map_chr(map, 1)),
#     bps = as.numeric(purrr::map_chr(map, 2)))

# }

# kraken %<>%
#   dplyr::mutate(
#     map_tib = stringr::str_split(map, " "),
#     map_tib = furrr::future_map(map_tib, map_to_tibble))

get_taxid <- function(taxa, kraken_taxids) {

  last_taxa <- stringr::str_split(taxa, "\\|")[[1]]
  last_taxa <- last_taxa[length(last_taxa)]
  taxa_word <- stringr::str_remove(last_taxa, regex("^[d|p|c|o|f|g|s]__"))
  out <- kraken_taxids %>%
    dplyr::filter(only_taxa == taxa_word)
  if (nrow(out) == 0) NA_real_ else min(out$taxid)

}

ksummary %<>%
  dplyr::mutate(
    taxid = purrr::map_dbl(taxa, get_taxid, kraken_taxids)) %>%
  dplyr::select(starts_with("tax"), everything())

check_hits_blast <- function(asv, taxa, taxid,
  blast, ksummary, kraken_taxids, min_prob = .1, max_dist = 2) {

  print(asv)

  if (is.null(blast)) {
    blast <- tibble::tibble()
  } else {
    blast %<>%
      dplyr::mutate(
        taxa = stringr::str_split(stitle, " "),
        taxa = purrr::map_chr(taxa, 1))

    nmatches <- nrow(blast)
    blast %<>%
      dplyr::count(taxa) %>%
      dplyr::mutate(prob = n / nmatches) %>%
      dplyr::filter(prob >= min_prob)
  }

  if (nrow(blast) > 0) {
    # get kraken id
    tid <- taxid
    ksummary_id <- ksummary %>%
      dplyr::filter(taxid == tid)

    # get blast id based on the genus
    ksummary_genus <- ksummary %>%
      dplyr::filter(!is.na(genus))
    ksummary_genus2 <- stringr::str_remove(ksummary_genus$genus, "g__")

    mmat <- stringdist::stringdistmatrix(blast$taxa,
      ksummary_genus2, method = "osa")
    rownames(mmat) <- blast$taxa

    mmat_search <- apply(mmat, 1, function(x)which(x < max_dist),
      simplify = FALSE)

    if (is.list(mmat_search)) {
      blast %<>%
        dplyr::mutate(
          taxonomy = mmat_search,
          taxonomy = purrr::map(taxonomy, ~ ksummary_genus[., ]),
          taxonomy = purrr::map(taxonomy, dplyr::select, -taxa, -species,
            -taxid),
          taxonomy = purrr::map(taxonomy, dplyr::distinct))
    } else {

      blast %<>%
        dplyr::mutate(
          taxonomy = list(ksummary_genus[mmat_search, ]),
          taxonomy = purrr::map(taxonomy, dplyr::select, -taxa, -species,
            -taxid),
          taxonomy = purrr::map(taxonomy, dplyr::distinct))

    }

    blast %<>%
      dplyr::mutate(
        prob = prob / sum(prob),
        d_hit = purrr::map_lgl(taxonomy,
          ~ ifelse(nrow(.) == 0, FALSE, .$domain == ksummary_id$domain)),
        p_hit = purrr::map_lgl(taxonomy,
          ~ ifelse(nrow(.) == 0, FALSE, .$phylum == ksummary_id$phylum)),
        c_hit = purrr::map_lgl(taxonomy,
          ~ ifelse(nrow(.) == 0, FALSE, .$class == ksummary_id$class)),
        o_hit = purrr::map_lgl(taxonomy,
          ~ ifelse(nrow(.) == 0, FALSE, .$order == ksummary_id$order)),
        f_hit = purrr::map_lgl(taxonomy,
          ~ ifelse(nrow(.) == 0, FALSE, .$family == ksummary_id$family)),
        g_hit = purrr::map_lgl(taxonomy,
          ~ ifelse(nrow(.) == 0, FALSE, .$genus == ksummary_id$genus)),
        across(
          where(is.logical),
          list(~ ifelse(is.na(.), 0, as.numeric(.))), .names = "{.col}"))

    blast %>%
      dplyr::select(ends_with("hit")) %>%
      colMeans()

  } else {

    rep(0, 6) %>%
      rlang::set_names(c("d_hit", "p_hit", "c_hit", "o_hit", "f_hit", "g_hit"))

  }

}


blast %<>%
  tidyr::nest(blast_pred = -c(qseqid)) %>%
  dplyr::rename(asv = qseqid)

kraken %<>%
  dplyr::left_join(blast, by = "asv")

kraken %<>%
  dplyr::mutate(
    hit_vec = purrr::pmap(list(asv, taxa, taxid, blast_pred),
      check_hits_blast, ksummary, kraken_taxids))

# outputs

# - accuracy hits
hits <- bind_rows(kraken$hit_vec) %>%
  dplyr::mutate(
    asv = kraken$asv,
    taxa = kraken$taxa,
    taxid = kraken$taxid) %>%
  dplyr::select(asv, taxid, tidyselect::everything())

# - rowData - asv and ksummary
rdata <- kraken %>%
  dplyr::select(asv, taxid) %>%
  dplyr::inner_join(ksummary, by = "taxid") %>%
  dplyr::select(asv, taxa, taxid, tidyselect::everything())


fs::dir_create(dir(arguments$taxa_table))
rdata %>%
  qs::qsave(arguments$taxa_table)

fs::dir_create(dir(arguments$hits_table))
hits %>%
  qs::qsave(arguments$hits_table)

