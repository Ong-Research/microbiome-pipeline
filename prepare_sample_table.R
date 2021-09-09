#!/ua/rwelch/miniconda3/envs/wiscdust/bin/Rscript

#' Creates a sample table with columns `batch|file_key|R1_file|R2_file` based on
#' the contents of a predefined directory
#' By default, it will look for the contents of the `data` directory and save
#' the table into `./samples.tsv`
#' @author rwelch

"Prepares a table to process the sequence files into an ASV table

Usage:
  prepare_sample_table.R [--outfile=<outfile> --input_dir=<inputdir>]
  prepare_sample_table.R (-h|--help)
  prepare_sample_table.R --version" -> doc
  
library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt(doc, args = my_args, version = "prepare sample table v1")

if (is.null(arguments$input_dir)) arguments$input_dir <- "data"
if (is.null(arguments$outfile)) arguments$outfile <- "samples.tsv"

stopifnot(dir.exists(arguments$input_dir))

message("generating sample table from ", arguments$input_dir)

library(magrittr, quietly = TRUE)
library(tidyverse, quietly = TRUE)

all_files <- list.files(arguments$input_dir, pattern = "fastq",
  recursive = TRUE)
  
# only keep fastq.gz files
all_files %<>%
  stringr::str_subset(stringr::regex("fastq.gz$"))

out <- tibble::tibble(filename = all_files) %>%
  dplyr::mutate(
    batch = stringr::str_split(filename, "\\/"),
    batch = purrr::map_chr(batch, 1),
    end = dplyr::if_else(stringr::str_detect(filename, "R1"), "end1", "end2"),
    key = stringr::str_split(filename, stringr::regex("_R[1|2]")),
    key = purrr::map_chr(key, 1)) %>%
  tidyr::pivot_wider(names_from = end, values_from = filename,
    values_fill = NA_character_) %>%
  na.omit() %>%
  dplyr::mutate(
    key = glue::glue("sample_{id}", id = dplyr::row_number()),
    end1 = file.path(arguments$input_dir, end1),
    end2 = file.path(arguments$input_dir, end2))

print(dplyr::count(out, batch))

message("saving sample table in ", arguments$outfile)

out %>%
  readr::write_tsv(arguments$outfile)