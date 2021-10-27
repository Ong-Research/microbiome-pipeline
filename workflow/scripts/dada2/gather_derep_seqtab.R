#!/usr/local/bin/Rscript

#' `dada2::makeSequenceTable` wrap
#' @param asv_file name of the ASV table file
#' @param summary_file name of the summary file with the number of reads per
#'   step
#' @param derep_file names of the file(s) with the dereplicated and merged # of
#'   reads vectors
#' @author rwelch2

"Gather ASV table

Usage:
gather_derep_seqtab.R [<asv_file> <summary_file>] [<derep_file> ...] [--batch=<batch> --log=<logfile> --config=<cfile>]
gather_derep_seqtab.R (-h|--help)
gather_derep_seqtab.R --version

Options:
-h --help    show this screen
--log=<logfile>    name of the log file [default: gather_derep_seqtab.log]
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "gather dereplicated sequence table V1")

if (!interactive()) {
  fs::dir_create(dirname(arguments$log))
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$batch <- "batch2018"
  arguments$asv_file <- "batch2018_asv.qs"
  arguments$derep_file <- list.files(file.path("output", "dada2",
    "merge", arguments$batch), full.names = TRUE)

}

message("arguments")
print(arguments)

message("info")
info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(dada2)
library(qs)
library(yaml)

derep_files <- arguments$derep_file

stopifnot(any(file.exists(derep_files)), file.exists(arguments$config))

# Identify and remove empty files.
# These are placeholder files that were created
# to appease snakemake after all reads were
# filtered out from the input file back in filter_and_trim.
# stop if we have no files left after checking for empty.
derep_tb <- tibble(derep_file = derep_files) %>%
  dplyr::mutate(
    size = purrr::map_dbl(derep_file, ~ file.info(.x)$size),
    key = purrr::map(derep_file, basename),
    key = stringr::str_remove(key, "_asv.qs"))

count_nonempty <- derep_tb %>%
  dplyr::filter(size > 0) %>%
  nrow()

stopifnot(count_nonempty > 0)

# read ASV files
read_merger <- function(filename, size) {
  if (size == 0) {
    NULL
  } else {
    qs::qread(filename)
  }
}

# changed the commands to be all inside the same mutate to improve readability
derep_tb %<>%
  dplyr::mutate(
    derep_merger = purrr::map2(derep_file, size, read_merger),
    mergers = purrr::map(derep_merger, "merge"),
    dada_fwd = purrr::map(derep_merger, "dada_fwd"),
    dada_bwd = purrr::map(derep_merger, "dada_bwd"))


config <- yaml::read_yaml(arguments$config)

stopifnot(file.exists(config$sample_table))

sample_tb <- readr::read_tsv(config$sample_table) %>%
  dplyr::left_join(derep_tb)

# now do some checking to make sure we only have
# the correct batch in the input
check_batches <- sample_tb %>%
  dplyr::filter(!is.na(derep_file)) %>%
  dplyr::count(batch)

if (nrow(check_batches) > 1) {
  message("User supplied input derep files from multiple batch(es):")
  message(paste(sprintf("%s: %d", check_batches$batch, check_batches$n),
    collapse = "\n"))
  message("We will only look at the requested batch: ", arguments$batch)
}

# remove other batches, but keep empty files for summary
sample_tb %<>%
  dplyr::filter(batch == arguments$batch)
nonempty_tb <- sample_tb %>%
  dplyr::filter(size > 0)
mergers <- pluck(nonempty_tb, "mergers") %>%
  setNames(pluck(nonempty_tb, "key"))
dada_fwd <- pluck(nonempty_tb, "dada_fwd")
dada_fwd <- pluck(nonempty_tb, "dada_bwd")

message("creating sequence table")
seqtab <- dada2::makeSequenceTable(mergers)

qs::qsave(seqtab, arguments$asv_file)

message("summarizing results")

## get N reads
## include lost samples with 0s
get_nreads <- function(x) {

  if (!is.null(x)) {
    sum(dada2::getUniques(x))
  } else {
    NA
  }
}

track <- sample_tb %>%
  dplyr::mutate(
    denoised = purrr::map(dada_fwd, get_nreads),
    merged = purrr::map(mergers, get_nreads)) %>%
  tidyr::unnest(c(denoised, merged)) %>%
  dplyr::transmute(samples = key, denoised, merged) %>%
  dplyr::mutate(across(c(denoised, merged), ~replace_na(.x, 0)))

track %>%
  readr::write_tsv(arguments$summary_file)

message("Done! summary file at ", arguments$summary_file)
