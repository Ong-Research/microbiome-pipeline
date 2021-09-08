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
gather_derep_seqtab.R [<asv_file> <summary_file>] [<derep_file> ...] [--log=<logfile> --config=<cfile>]
gather_derep_seqtab.R (-h|--help)
gather_derep_seqtab.R --version

Options:
-h --help    show this screen
--log=<logfile>    name of the log file [default: filter_and_trim.log]
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "gather dereplicated sequence table V1")


log_file <- file(arguments$log, open = "wt")
sink(log_file, type = "message")

print(arguments)
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

config <- yaml::read_yaml(arguments$config)

stopifnot(file.exists(config$sampletable))
sample_names <- readr::read_tsv(config$sampletable)[[1]]

derep_mergers <- purrr::map(derep_files, qs::qread)

mergers <- purrr::map(derep_mergers, "merge")
names(mergers) <- sample_names

dada_fwd <- purrr::map(derep_mergers, "dada_fwd")
dada_bwd <- purrr::map(derep_mergers, "dada_bwd")

message("creating sequence table")
seqtab <- dada2::makeSequenceTable(mergers)

qs::qsave(seqtab, arguments$asv_file)

message("summarizing results")

## get N reads
get_nreads <- function(x) sum(dada2::getUniques(x))

track <- cbind(
    sapply(dada_fwd, get_nreads), sapply(mergers, get_nreads))
colnames(track) <- c("denoised", "merged")

track %>%
  as.data.frame() %>%
  tibble::as_tibble(rownames = "samples") %>%
  readr::write_tsv(arguments$summary_file)

message("Done! summary file at ", arguments$summary_file)

close(log_file)
