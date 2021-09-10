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
filter_and_trim.R [<filter_end1> <filter_end2> <summary_file>] [<sample_name> --end1=<end1> --end2=<end2>] [--batch=<batch> --log=<logfile> --config=<cfile>]
filter_and_trim.R (-h|--help)
filter_and_trim.R --version

Options:
-h --help    show this screen
--end1=<end1>    name of the R1 end fastq.gz file
--end2=<end2>    name of the R2 end fastq.gz file
--log=<logfile>    name of the log file [default: logs/filter_and_trim.log]
--batch=<batch>    name of the batch if any to get the filter and trim parameters
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args, version = "filter_and_trim V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$batch <- "dust_dec2018"
  arguments$sample_name <- "sample_20"
  arguments$end1 <-
    "data/dust_dec2018/190114_C75PR/203_S19_L001_R1_001.fastq.gz"
  arguments$end2 <-
    "data/dust_dec2018/190114_C75PR/203_S19_L001_R2_001.fastq.gz"

  arguments$filter_end1 <-
    "filtered_L001_R1_001.fastq.gz"
  arguments$filter_end2 <-
    "filtered_L001_R2_001.fastq.gz"
}

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
config <- yaml::read_yaml(arguments$config)
config <- config[["filter_and_trim"]]
print(config)

if (!is.null(arguments$batch)) {
  stopifnot(arguments$batch %in% names(config))
  config <- config[[arguments$batch]]
} else {
  nms <- c("truncQ", "truncLen", "trimLeft", "trimRight",
    "maxLen", "minLen", "maxN", "minQ", "maxEE")
  if (any(names(config) %in% nms)) {
    warning("will use first element instead")
    config <- config[[1]]
  }

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
  dplyr::mutate(
    end1 = arguments$end1,
    end2 = arguments$end2) %>%
  readr::write_tsv(arguments$summary_file)
