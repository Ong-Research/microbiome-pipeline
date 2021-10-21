#!/usr/local/bin/Rscript

#' Subtract p% of the negative controls reads of the ASV matrix
#' @author rwelch2

"Filters the low quality ASVs

Usage:
filter_asvs.R [<asv_matrix_qc> <seqlength_fig> <abundance_fig>] [<asv_matrix_file> <negcontrol_file>] [--lysis --log=<logfile> --config=<cfile> --cores=<cores>]
filter_asvs.R (-h|--help)
filter_asvs.R --version

Options:
-h --help    show this screen
--log=<logfile>    name of the log file [default: filter_and_trim.log]
--config=<cfile>    name of the yaml file with the parameters [default: ./config/config.yaml]
--cores=<cores>    number of CPUs for parallel processing [default: 24]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "filter ASVs V1")

if (!interactive()) {
  fs::dir_create(dirname(arguments$log))
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$asv_matrix_qc <- "asv.qs"
  arguments$seqlength_fig <- "sl.png"
  arguments$abundance_fig <- "sa.png"
  arguments$asv_matrix_file <- "output/dada2/remove_chim/asv_mat_wo_chim.qs"
  arguments$negcontrol_file <- "data/negcontrol.qs"

}


message("arguments")
print(arguments)

info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))


message("loading packages")
library(magrittr)
library(tidyverse)
library(dada2)
library(qs)
library(yaml)

stopifnot(
  file.exists(arguments$asv_matrix_file),
  file.exists(arguments$config),
  file.exists(arguments$negcontrol_file))

message("starting with asvs in ", arguments$asv_matrix_file)

seqtab <- qs::qread(arguments$asv_matrix_file)
config <- yaml::read_yaml(arguments$config)[["qc_parameters"]]

message("filtering samples by negative controls")
message("using neg. control file ", arguments$negcontrol_file)
message("removing counts >= ", config[["negctrl_prop"]],
  " sum(neg_controls)")

neg_controls <- qs::qread(arguments$negcontrol_file)

if (arguments$lysis) {

  neg_controls %<>%
    dplyr::mutate(
      kits = map2(kits, lysis, c))

}

neg_controls %<>%
  dplyr::select(batch, key, kits)

# only look at samples that
# are in the input matrix
neg_controls %<>%
  dplyr::filter(key %in% rownames(seqtab))

subtract_neg_control <- function(name, neg_controls, seqtab, prop) {

  negs <- neg_controls

  negs <- negs[negs %in% rownames(seqtab)]
  out_vec <- seqtab[name, ]

  if (length(negs) > 1) {

    negs <- seqtab[negs, ]
    neg_vec <- colSums(negs)
    out_vec <- out_vec - prop * neg_vec

  } else if (length(negs) == 1) {

    neg_vec <- seqtab[negs, ]
    out_vec <- out_vec - prop * neg_vec
  }

  if (any(out_vec < 0)) {
    out_vec[out_vec < 0] <- 0
  }
  return(out_vec)

}

neg_controls %<>%
  dplyr::mutate(sample_vec = purrr::map2(key, kits,
    safely(subtract_neg_control),
    seqtab, config[["negctrl_prop"]]))

to_remove <- dplyr::filter(neg_controls,
  purrr::map_lgl(sample_vec, ~ !is.null(.$error)))

if (nrow(to_remove) > 0) {
  message(sprintf("removed %d samples:\n", nrow(to_remove)),
    stringr::str_c(to_remove$key, collapse = ", "))
}

neg_controls %<>%
  dplyr::filter(purrr::map_lgl(sample_vec, ~ is.null(.$error))) %>%
  dplyr::mutate(
    sample_vec = purrr::map(sample_vec, "result"))

sample_names <- neg_controls %>%
  dplyr::pull(key)

seqtab_samples <- purrr::reduce(
  dplyr::pull(neg_controls, sample_vec), dplyr::bind_rows)
seqtab_samples %<>%
  as.matrix() %>%
  set_rownames(sample_names)

seqtab_nc <- seqtab[! rownames(seqtab) %in% sample_names, ]
seqtab_new <- floor(rbind(seqtab_samples, seqtab_nc))

# Length of sequences
message("filtering sequences by length")

seq_lengths <- nchar(dada2::getSequences(seqtab_new))

l_hist <- as.data.frame(table(seq_lengths)) %>%
  tibble::as_tibble() %>%
  rlang::set_names("length", "freq")

l_hist <- l_hist %>%
  ggplot(aes(x = length, y = freq)) +
    geom_col() +
    labs(title = "Sequence Lengths by SEQ Count") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, hjust = 1,
        vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10))

ggsave(
  filename = arguments$seqlength_fig,
  plot = l_hist,
  width = 8, height = 4, units = "in")


table2 <- tapply(colSums(seqtab_new), seq_lengths, sum)
table2 <- tibble::tibble(
  seq_length = names(table2),
  abundance = table2)

most_common_length <- dplyr::top_n(table2, 1, abundance) %>%
  dplyr::pull(seq_length) %>%
  as.numeric()

table2 <- table2 %>%
  ggplot(aes(x = seq_length, y = log1p(abundance))) +
    geom_col() +
    labs(title = "Sequence Lengths by SEQ Abundance") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, hjust = 1,
        vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10))

ggsave(
  filename = arguments$abundance_fig,
  plot = table2,
  width = 8, height = 4, units = "in")


max_diff <- config[["max_length_variation"]]

message("most common length: ", most_common_length)
message("removing sequences outside range ",
  " < ", most_common_length - max_diff, " or > ",
  most_common_length + max_diff)

right_length <- abs(seq_lengths - most_common_length) <= max_diff
seqtab_new <- seqtab_new[, right_length]


total_abundance <- sum(colSums(seqtab_new))
min_reads_per_asv <- ceiling(config[["low_abundance_perc"]] / 100 *
  total_abundance)

seqtab_abundance <- colSums(seqtab_new)

message("removing ASV with < ", min_reads_per_asv, " reads")
message("in total ", sum(seqtab_abundance < min_reads_per_asv))
seqtab_new <- seqtab_new[, seqtab_abundance >= min_reads_per_asv]


message("saving files in ", arguments$asv_matrix_qc)
qs::qsave(seqtab_new, arguments$asv_matrix_qc)
