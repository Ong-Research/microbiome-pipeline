#!/usr/local/bin/Rscript

#' Clean the metadata file into a tsv file
#' @author rwelch2

"Clean metadata into a tsv file

Usage:
clean_metadata.R [<tsv_file>] [--asv=<asv_file> --meta=<meta_file>] [--log=<logfile>]
clean_metadata.R (-h|--help)
clean_metadata.R --version

Options:
-h --help    show this screen
--asv=<asv_file>    ASV matrix file
--meta=<meta_file>    Metadata file
--log=<logfile>    name of the log file [default: logs/filter_and_trim.log]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "clean metadata V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {
  arguments$asv <- "output/dada2/after_qc/asv_mat_wo_chim.qs"
  arguments$meta <- "output/predada2/meta.qs"
}

print(arguments)

info <- Sys.info();

message("loading packages")
library(magrittr)
library(tidyverse)
library(qs)

stopifnot(
  file.exists(arguments$asv),
  file.exists(arguments$meta)
)

asv <- qs::qread(arguments$asv)
meta <- qs::qread(arguments$meta)

# need to fix parse_taxa to return a tibble of length 
# equal to all the sequences and not only the sequences with known information

icra_data <- meta$data[[3]] %>%
  dplyr::filter(months <= 24) %>%
  dplyr::select(
    tidyselect::ends_with("sid"),
    months,
    icra_3_runny_nose,
    icra_7_cough,
    icra_13_wheezing,
    icra_28_eczema) %>%
  tidyr::pivot_wider(ends_with("sid"), names_from = months,
    values_from = starts_with("icra"))
  
pre_data <- meta$data[[4]] %>%
  dplyr::select(-pilot, -farm)

cdata <- tibble::tibble(key = rownames(asv)) %>%
  dplyr::left_join(meta$data[[1]], by = "key") %>%
  dplyr::left_join(meta$data[[2]], by = "baby_sid") %>%
  dplyr::left_join(pre_data, by = c("baby_sid", "mom_sid")) %>%
  dplyr::left_join(icra_data, by = c("baby_sid", "mom_sid"))


fs::dir_create(dirname(arguments$tsv_file))
readr::write_tsv(cdata, arguments$tsv_file)
