#!/usr/local/bin/Rscript

#' `dada2::mergeSequenceTables` and `dada2::removeBimeraDenovo` wrapper
#' @author rwelch2

"Merge sequence tables and remove chimeras

Usage:
remove_chimeras.R [<asv_merged_file> <summary_file> <fig_file>] [<asv_file> ...] [--log=<logfile> --config=<cfile> --cores=<cores>]
remove_chimeras.R (-h|--help)
remove_chimeras.R --version

Options:
-h --help    show this screen
--log=<logfile>    name of the log file [default: remove_chimeras.log]
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]
--cores=<cores>    number of CPUs for parallel processing [default: 24]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "remove chimeras V1")

if (!interactive()) {
  fs::dir_create(dirname(arguments$log))
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$fig_file <- "matspy.png"
  arguments$asv_merged_file <- "asv.qs"
  arguments$summary_file <- "summary.qs"
  arguments$asv_file <- list.files(file.path("output", "dada2",
    "asv_batch"), full.names = TRUE, pattern = "asv")

}


message("arguments")
print(arguments)

## merges different ASV tables, and removes chimeras

message("info")
info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(dada2)
library(qs)
library(yaml)
library(ComplexHeatmap)

stopifnot(any(file.exists(arguments$asv_file)))

if (!is.null(arguments$config)) stopifnot(file.exists(arguments$config))

config <- yaml::read_yaml(arguments$config)
sample_table <- readr::read_tsv(config$sample)
config <- config$remove_chimeras

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

fs::dir_create(dirname(arguments$fig_file))

outmat <- seqtab[, seq_len(min(1e4, ncol(seqtab)))]
outmat <- ifelse(outmat > 0, 1, 0)
cols <- structure(c("black", "white"), names = c("1", "0"))


# put row annotation in same order as matrix
# will remove samples that were completely empty
# after filter_and_trim
anno_df <- sample_table %>%
    select(batch, key) %>%
    as.data.frame() %>%
    column_to_rownames("key")
anno_df <- anno_df[rownames(outmat), , drop = F]

annot <- ComplexHeatmap::rowAnnotation(
  df = anno_df,
  annotation_legend_param = list(
    batch = list(direction = "horizontal")))

png(filename = arguments$fig_file, width = 10, height = 8, units = "in",
  res = 1200)
hm <- ComplexHeatmap::Heatmap(outmat,
  left_annotation = annot,
  col = cols, name = "a", show_row_dend = FALSE, show_row_names = FALSE,
  show_column_dend = FALSE, show_column_names = FALSE,
  cluster_columns = FALSE, cluster_rows = FALSE,
  show_heatmap_legend = FALSE, use_raster = TRUE,
  column_title = "Top 10K ASVs",
  heatmap_legend_param = list(direction = "horizontal"))
draw(hm, annotation_legend_side = "top", heatmap_legend_side = "top",
  merge_legend = TRUE)
dev.off()
