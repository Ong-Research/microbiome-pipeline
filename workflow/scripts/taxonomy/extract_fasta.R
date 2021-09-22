#!/usr/local/bin/Rscript

#' Extracts fasta sequences from the ASV matrices

"Extract fasta sequences

Usage:
extract_fasta.R [<fasta_file>] [<asv_mat_file>] [--prefix=<pref> --log=<logfile>]
extract_fasta.R (-h|--help)
extract_fasta.R --version

Options:
--prefix=<pref>    prefix to name the sequences [default: asv]
--log=<logfile>    name of the log file [default: extract_fa.log]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "extract fasta V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$fasta_file <- "out.fasta"
  arguments$asv_mat_file <-
    "output/dada2/after_qc/asv_mat_wo_chim.qs"

}
   
info <- Sys.info();

print(arguments)
print(stringr::str_c(names(info), " : ", info, "\n"))

stopifnot(file.exists(arguments$asv_mat_file))

message("loading packages")
library(magrittr)
library(tidyverse)
library(Biostrings)
library(qs)

asvs <- qs::qread(arguments$asv_mat_file)
sequences <- colnames(asvs)
names(sequences) <- stringr::str_c(arguments$prefix,
  seq_along(sequences), sep = "_")

fs::dir_create(dirname(arguments$fasta_file))

sequences <- Biostrings::DNAStringSet(sequences)
Biostrings::writeXStringSet(sequences, filepath = arguments$fasta_file)
