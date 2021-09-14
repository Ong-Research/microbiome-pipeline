#!/usr/local/bin/Rscript

#' Collect multiple tsv files into one provides they have the same column names

"Collect tsv files into one

Usage:
collect_tsv_files.R [<outfile>] [<infiles> ...] [--log=<logfile>]
collect_tsv_files.R (-h|--help)
collect_tsv_files.R --version

Options:
--log=<logfile>    name of the log file [default: collect.log]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "collect tsv files V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$outfile <- "out.tsv"
  arguments$infiles <- list.files("output/dada2/summary",
    full.names = TRUE)

}
   
info <- Sys.info();

print(arguments)
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")

library(magrittr)
library(tidyverse)
library(vroom)

infiles <- purrr::map(arguments$infiles, vroom::vroom)

out <- bind_rows(infiles)

fs::dir_create(dirname(arguments$outfile))
readr::write_tsv(out, file = arguments$outfile)
