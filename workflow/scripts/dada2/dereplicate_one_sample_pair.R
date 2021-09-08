#!/usr/local/bin/Rscript

#' Wrapper around `dada2::dada2::derepFastq`, `dada2::dada` and
#' `dada2::mergePairs` functions
#' For details on the meaning of the parameters use
#' R -e '?dada2::derepFastq'
#' R -e '?dada2::dada'
#' R -e '?dada2::mergePairs'
#' @param sample_merge_file name of the output merged ASV vector
#' @param sample_name name of the sample
#' @param end1_file name of the filtered end1 file
#' @param end2_file name of the filtered end2 file
#' @author rwelch

#' This script is a bit more elaborated, because instead of wrapping only one
#' dada2 method, this wraps the following:
#'
#' * dada2::derepFastq
#' * dada2::dada
#' * dada2::mergePairs

"Dereplicate one sample pair

Usage:
dereplicate_one_sample_pair.R  [<sample_merge_file>] [<sample_name> <end1_file> <end2_file> --end1_err=<end1_err> --end2_err=<end2_err>] [--log=<logfile> --batch=<batch> --config=<cfile>]
dereplicate_one_sample_pair.R (-h|--help)
dereplicate_one_sample_pair.R --version

Options:
--end1_err=<end1>    name of the R1 error rates matrix
--end2_err=<end2>    name of the R2 error rates matrix
--log=<logfile>    name of the log file [default: dereplicate.log]
--batch=<batch>    name of the batch if any to get the filter and trim parameters
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "dereplicate one sample V1")

stopifnot(
  file.exists(arguments$end1_file),
  file.exists(arguments$end2_file),
  file.exists(arguments$end1_err),
  file.exists(arguments$end2_err))

if (!is.null(arguments$config)) stopifnot(file.exists(arguments$config))
  
log_file <- file(arguments$log, open = "wt")
sink(log_file)

info <- Sys.info();

print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")

library(magrittr)
library(dada2)
library(qs)
library(yaml)

filter_fwd <- arguments$end1_file
filter_bwd <- arguments$end2_file

message("loading error rates")
err_fwd <- qs::qread(arguments$end1_err)
err_bwd <- qs::qread(arguments$end2_err)

message("dereplicating filtered files")

dereplicate_merge <- function(filtfwd, filtbwd, err_fwd, err_bwd,
  min_overlap = 12, max_mismatch = 0) {

  message("de-replicating end1 file")
  derep_fwd <- dada2::derepFastq(filtfwd)
  dd_fwd <- dada2::dada(derep_fwd, err = err_fwd)

  message("de-replicating end2 file")
  derep_bwd <- dada2::derepFastq(filtbwd)
  dd_bwd <- dada2::dada(derep_bwd, err = err_bwd)

  message("merging pairs")
  merge <- dada2::mergePairs(dd_fwd, derep_fwd, dd_bwd, derep_bwd,
      minOverlap = min_overlap, maxMismatch = max_mismatch)

  out <- list(
    "dada_fwd" = dd_fwd,
    "dada_bwd" = dd_bwd,
    "merge" = merge)

  return(out)
}

config <- yaml::read_yaml(arguments$config)$merge_pairs

if (!is.null(arguments$batch)) {
  stopifnot(arguments$batch %in% names(config))
  config <- config[[arguments$batch]]
}


merge <- dereplicate_merge(filter_fwd, filter_bwd, err_fwd, err_bwd,
  min_overlap = as.numeric(config[["minOverlap"]]),
  max_mismatch = as.numeric(config[["maxMismatch"]]))

qs::qsave(merge, arguments$sample_merge_file)

close(log_file)
