# this script filters the ASV matrix after removing chimeras
#
# it does the following stuff:
# - 
source("renv/activate.R")
info <- Sys.info();

sink(snakemake@log[[1]])

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(dada2)
library(qs)

message("starting with asvs in ", snakemake@input[["seqtab"]])

seqtab <- qs::qread(snakemake@input[["seqtab"]])

message("filtering samples by negative controls")
message("using neg. control file ", snakemake@input[["negcontrol"]])
message("removing counts >= ", snakemake@config[["negctrl_prop"]],
  " sum(neg_controls)")
neg_controls <- readr::read_tsv(snakemake@input[["negcontrol"]])

neg_controls %<>%
  tidyr::nest(negs = c(neg_control))

subtract_neg_control <- function(name, neg_controls, seqtab, prop) {

  negs <- neg_controls %>%
    dplyr::pull(neg_control)

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
  dplyr::mutate(sample_vec = purrr::map2(name, negs,
    safely(subtract_neg_control),
    seqtab, snakemake@config[["negctrl_prop"]]))

to_remove <- dplyr::filter(neg_controls,
  purrr::map_lgl(sample_vec, ~ !is.null(.$error)))

if (nrow(to_remove) > 0) {
  message("removed the samples:\n",
    stringr::str_c(to_remove$name, collapse = "\n"))
}

neg_controls %<>%
  dplyr::filter(purrr::map_lgl(sample_vec, ~ is.null(.$error))) %>%
  dplyr::mutate(
    sample_vec = purrr::map(sample_vec, "result"))

sample_names <- neg_controls %>%
  dplyr::pull(name)

seqtab_samples <- purrr::reduce(
  dplyr::pull(neg_controls, sample_vec), dplyr::bind_rows)
seqtab_samples %<>%
  as.matrix() %>%
  set_rownames(sample_names)

seqtab_nc <- seqtab[! rownames(seqtab) %in% sample_names, ]
seqtab_new <- rbind(seqtab_samples, seqtab_nc)

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
      axis.text.x = element_text(angle = 90, hjust = 1,
        vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10))

ggsave(
  filename = snakemake@output[["plot_seqlength"]],
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
      axis.text.x = element_text(angle = 90, hjust = 1,
        vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10))

ggsave(
  filename = snakemake@output[["plot_seqabundance"]],
  plot = table2,
  width = 8, height = 4, units = "in")


max_diff <- snakemake@config[['max_length_variation']]

message("most common length: ", most_common_length)
message("removing sequences outside range ",
  " < ", most_common_length - max_diff, " or > ",
  most_common_length + max_diff)

right_length <- abs(seq_lengths - most_common_length) <= max_diff
seqtab_new <- seqtab_new[, right_length]


total_abundance <- sum(colSums(seqtab_new))
min_reads_per_asv <- ceiling(snakemake@config[["low_abundance_perc"]] / 100 *
  total_abundance)

seqtab_abundance <- colSums(seqtab_new)

message("removing ASV with < ", min_reads_per_asv, " reads")
message("in total ", sum(seqtab_abundance < min_reads_per_asv))
seqtab_new <- seqtab_new[, seqtab_abundance >= min_reads_per_asv]


message("saving files in ", snakemake@output[["seqtab_filt"]])
qs::qsave(seqtab_new, snakemake@output[["seqtab_filt"]])
