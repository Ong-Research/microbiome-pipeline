#!/usr/local/bin/Rscript

#' Summarize number of reads per processing step
#' @author rwelch2
#' 
"Summarized the # of reads per processing step

Usage:
summarize_nreads.R [<nreads_file> <nreads_fig> <preads_fig>] [<filt_summary_file> <derep_summary_file> <final_asv_mat>] [--log=<logfile>]
summarize_nreads.R (-h|--help)
summarize_nreads.R --version

Options:
-h --help    show this screen
--log=<logfile>    name of the log file [default: summarize_nreads.log]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "summarize # reads V1")

if (!interactive()) {
  fs::dir_create(dirname(arguments$log))
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$filt_summary_file <-
    "output/dada2/filtered/all_sample_summary.tsv"
  arguments$derep_summary_file <-
    "output/dada2/asv_batch/all_sample_summary.tsv"
  arguments$final_asv_mat <-
    "output/dada2/after_qc/asv_mat_wo_chim.qs"

  arguments$nreads_file <-
    "output/dada2/stats/Nreads_dada2.txt"
  arguments$nreads_fig <-
    "workflow/report/dada2qc/dada2steps_vs_abundance.png"
  arguments$preads_fig <-
    "workflow/report/dada2qc/dada2steps_vs_relabundance.png"

}

# Summarize numbers of reads per step, makes some plots

info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(vroom)
library(qs)

stats <- list()
stats[[1]] <- readr::read_tsv(arguments$filt_summary_file) %>%
  select(-end1, -end2)
stats[[2]] <- readr::read_tsv(arguments$derep_summary_file)
stats[[3]] <- qs::qread(arguments$final_asv_mat) %>%
  rowSums() %>%
  tibble::tibble(samples = names(.), nreads = .)

#stats <- purrr::reduce(stats, purrr::partial(dplyr::inner_join, by = "samples"))

# change to full join so we see input files that were lost during any step
stats <- purrr::reduce(stats, purrr::partial(dplyr::full_join, by = "samples"))
# fill in 0s for final nreads for samples that were previously lost
stats %<>%
  dplyr::mutate(nreads = replace_na(nreads, 0))

stats %>%
  readr::write_tsv(arguments$nreads_file)

message("making figures")

order_steps <- stats %>%
  dplyr::select(-samples) %>%
  names()

rel_stats <- stats %>%
  dplyr::mutate(
    dplyr::across(
      -samples, list(~ . / raw), .names = "{.col}"))

# remove empties from relative change plot
# they don't get plotted anyway
rel_stats %<>%
  dplyr::filter(!is.nan(nreads))

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
  filename = arguments$nreads_fig,
  plot = make_plot(stats) + labs("# reads") +
   scale_y_continuous(labels = scales::comma_format(1)),
  width = 6,
  height = 4,
  units = "in")

ggsave(
  filename = arguments$preads_fig,
  plot = make_plot(rel_stats) + labs(y = "relative change") +
    scale_y_continuous(labels = scales::percent_format(1)),
  width = 6,
  height = 4,
  units = "in")
