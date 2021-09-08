#!/usr/local/bin/Rscript

#' Wrapper around `dada2::filterAndTrim` function
#' For details on the meaning of the parameters use
#' R -e '?dada2::filterAndTrim'
#' @param filter_end1 name of the output end1 file
#' @param filter_end2 name of the output end2 file
#' @param summary_file name of the file where the summary is saved
#' @param sample_name name of the sample
#' @author rwelch2

"Filter and trim

Usage:
filter_and_trim.R [<filter_end1> <filter_end2> <summary_file>] [<sample_name> --end1=<end1> --end2=<end2>] [--log=<logfile> --batch=<batch> --config=<cfile>]
filter_and_trim.R (-h|--help)
filter_and_trim.R --version

Options:
-h --help    show this screen
--end1=<end1>    name of the R1 end fastq.gz file
--end2=<end2>    name of the R2 end fastq.gz file
--log=<logfile>    name of the log file [default: filter_and_trim.log]
--batch=<batch>    name of the batch if any to get the filter and trim parameters
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args, version = "filter_and_trim V1")


log_file <- file(arguments$log, open = "wt")
sink(log_file, type = "message")

print(arguments)

info <- Sys.info();

message("loading dada2")
library(magrittr)
library(dada2)
library(qs)
library(yaml)
library(fs)

stopifnot(file.exists(arguments$config),
  file.exists(arguments$end1), file.exists(arguments$end2))

print(stringr::str_c(names(info), " : ", info, "\n"))
config <- yaml::read_yaml(arguments$config)$filter_and_trim
print(config)

if (!is.null(arguments$batch)) {
  stopifnot(arguments$batch %in% names(config))
  config <- config[[arguments$batch]]
}

fs::dir_create(unique(dirname(arguments$filter_end1)))
fs::dir_create(unique(dirname(arguments$filter_end2)))

track_filt <- dada2::filterAndTrim(
  arguments$end1, arguments$filter_end1,
  arguments$end2, arguments$filter_end2,
  truncQ =    as.numeric(config[["truncQ"]]),
  truncLen =  as.numeric(config[["truncLen"]]),
  trimLeft =  as.numeric(config[["trimLeft"]]),
  trimRight = as.numeric(config[["trimRight"]]),
  maxLen =    as.numeric(config[["maxLen"]]),
  minLen =    as.numeric(config[["minLen"]]),
  maxN =      as.numeric(config[["maxN"]]),
  minQ =      as.numeric(config[["minQ"]]),
  maxEE =     as.numeric(config[["maxEE"]]),
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE)

row.names(track_filt) <- arguments$sample_name
colnames(track_filt) <- c("raw", "filtered")

track_filt %>%
  as.data.frame() %>%
  tibble::as_tibble(rownames = "samples") %>%
  readr::write_tsv(arguments$summary_file)

close(log_file)
