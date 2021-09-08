# Summarize numbers of reads per step, makes some plots
source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(qs)

stats <- list()
stats[[1]] <- readr::read_tsv(snakemake@input[["nreads_filtered"]])
stats[[2]] <- readr::read_tsv(snakemake@input[["nreads_dereplicated"]])
stats[[3]] <- readr::read_tsv(snakemake@input[["nreads_chim_removed"]])

stats <- purrr::reduce(stats, purrr::partial(dplyr::inner_join, by = "samples"))

stats %>%
  readr::write_tsv(snakemake@output[["nreads"]])

message("making figures")

order_steps <- stats %>%
  dplyr::select(-samples) %>%
  names()

rel_stats <- stats %>%
  dplyr::mutate(
    dplyr::across(
      -samples, list( ~ . / raw), .names = "{.col}"))

make_plot <- function(stats, summary_fun = median, ...) {

  order_steps <- stats %>%
    dplyr::select(-samples) %>%
    names()

  stats %<>%
    tidyr::pivot_longer(-samples, names_to = "step", values_to = "val") %>%
    dplyr::mutate(step = factor(step, levels = order_steps))

  summary <- stats %>%
    dplyr::group_by(step) %>%
    dplyr::summarize(
      val = summary_fun(val, ...), .groups = "drop")

  stats %>%
    ggplot(aes(step, val)) + geom_boxplot() +
    geom_point(alpha = 0.25, shape = 21) +
    geom_line(aes(group = samples), alpha = 0.25) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.text.y = ggplot2::element_text(angle = -90, size = 10),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()) +
    geom_line(data = summary, aes(group = 1), linetype = 2, colour = "red") +
    labs(x = "step")

}

ggsave(
  filename = snakemake@output[["fig_step"]],
  plot = make_plot(stats) + labs("# reads"),
  width = 6,
  height = 4,
  units = "in")

ggsave(
  filename = snakemake@output[["fig_step_rel"]],
  plot = make_plot(rel_stats) + labs(y = "relative change"),
  width = 6,
  height = 4,
  units = "in")
