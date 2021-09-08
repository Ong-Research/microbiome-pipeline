#!/usr/local/bin/Rscript

#' `dada2::mergeSequenceTables` and `dada2::removeBimeraDenovo` wrapper
#' @author rwelch2

"Merge sequence tables and remove chimeras

Usage:
remove_chimeras.R [<asv_merged_file> <summary_file>] [<asv_file> ...] [--log=<logfile> --config=<cfile> --cores=<cores>]
remove_chimeras.R (-h|--help)
remove_chimeras.R --version

Options:
-h --help    show this screen
--log=<logfile>    name of the log file [default: filter_and_trim.log]
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]
--cores=<cores>    number of CPUs for parallel processing [default: 4]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "remove chimeras V1")


log_file <- file(arguments$log, open = "wt")
sink(log_file, type = "message")

print(arguments)

## merges different ASV tables, and removes chimeras

info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library
library(dada2)
library(qs)
library(yaml)

stopifnot(any(file.exists(arguments$asv_file)))

if (!is.null(arguments$config)) stopifnot(file.exists(arguments$config))

config <- yaml::read_yaml(arguments$config)$remove_chimeras

seqtab_list <- purrr::map(arguments$asv_file, qs::qread)

if (length(seqtab_list) > 1) {
  seqtab_all <- dada2::mergeSequenceTables(tables = seqtab_list)
} else {
  seqtab_all <- seqtab_list[[1]]
}



# Remove chimeras
message("removing chimeras")
seqtab <- dada2::removeBimeraDenovo(
  seqtab_all,
  method = config[["chimera_method"]],
  minSampleFraction = config[["minSampleFraction"]],
  ignoreNNegatives = config[["ignoreNNegatives"]],
  minFoldParentOverAbundance = config[["minFoldParentOverAbundance"]],
  allowOneOf = config[["allowOneOf"]],
  minOneOffParentDistance = config[["minOneOffParentDistance"]],
  maxShift = config[["maxShift"]],
  multithread = as.numeric(arguments$cores))

fs::dir_create(dirname(arguments$asv_merged_file))
qs::qsave(seqtab, arguments$asv_merged_file)

out <- tibble::tibble(samples = row.names(seqtab),
  nonchim = rowSums(seqtab))

fs::dir_create(dirname(arguments$summary_file))
out %>%
  readr::write_tsv(arguments$summary_file)

close(log_file)
