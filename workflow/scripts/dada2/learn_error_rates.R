#!/usr/local/bin/Rscript

#' Wrapper around `dada2::learnErrors` function
#' For details on the meaning of the parameters use
#' R -e '?dada2::learnErrors'
#' @param filtered name of filtered files
#' @author rwelch2

"Learn error rates

Usage:
learn_error_rates.R [--error_rates=<matrix_file> --plot_file=<pfile>] [<filtered> ...] [--log=<logfile> --batch=<batch> --config=<cfile> --cores=<cores>]
learn_error_rates.R (-h|--help)
learn_error_rates.R --version

Options:
-h --help    show this screen
--error_rates=<matrix_file>    name of the file with the learned error rates matrix
--plot_file=<pfile>    name of the file to save the diagnostic plot
--log=<logfile>    name of the log file [default: error_rates.log]
--batch=<batch>    name of the batch if any to get the filter and trim parameters
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]
--cores=<cores>    number of CPUs for parallel processing [default: 4]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args, version = "error_rates V1")

stopifnot(any(file.exists(arguments$filtered)))

if (!is.null(arguments$config)) {
  stopifnot(file.exists(arguments$config))
}

log_file <- file(arguments$log, open = "wt")

info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(dada2)
library(ggplot2)
library(qs)
library(yaml)
library(fs)

config <- yaml::read_yaml(arguments$config)$error_rates
print(config)

if (!is.null(arguments$batch)) {
  stopifnot(arguments$batch %in% names(config))
  config <- config[[arguments$batch]]
}

print("computing error rates")
errs <- dada2::learnErrors(
  arguments$filtered, nbases = as.numeric(config$learn_nbases),
  multithread = as.numeric(arguments$cores), randomize = TRUE)
  
fs::dir_create(dirname(arguments$error_rates))
qs::qsave(errs, file = arguments$error_rates)

print("plotting diagnostics")
err_plot <- dada2::plotErrors(errs, nominalQ = TRUE)

fs::dir_create(dirname(arguments$plot_file))
ggplot2::ggsave(
  filename = arguments$plot_file,
  plot = err_plot,
  width = 20,
  height = 20,
  units = "cm")

close(log_file)
