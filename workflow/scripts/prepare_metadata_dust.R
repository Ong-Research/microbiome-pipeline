#!/ua/rwelch/miniconda3/envs/wiscdust/bin/Rscript

#' Updates old metadata file to use the same keys as in the `./samples.tsv` file
#' @author rwelch

"Prepares a table to match negative controls to samples

Usage:
prepare_metadata_dust.R [<outfile>] [--samples=<sample_table> --xlsx_file=<xlsx_file>]
prepare_metadata_dust.R (-h|--help)
prepare_metadata_dust.R --version

Options:
--samples=<sample_table>   Name of the file with the samples table.
--xlsx_file=<xlsx_file>    Name of an xlsx table with metadata." -> doc
  
library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt(doc, args = my_args,
  version = "prepare metadata dust v1")

if (is.null(arguments$samples)) {
  arguments$samples <- "samples.tsv"
}
if (is.null(arguments$outfile)) arguments$outfile <- "data/meta.qs"

if (interactive()) {
  arguments$xlsx_file <-
    "data/wisc_dust_16S_samplelist.xlsx"
}


stopifnot(file.exists(arguments$samples))

library(magrittr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(readxl, quietly = TRUE)
library(qs, quietly = TRUE)


sample_table <- readr::read_tsv(arguments$samples)

xls_table <- readxl::read_xlsx(arguments$xlsx_file,
  sheet = readxl::excel_sheets(arguments$xlsx_file)[1],
  na = c("N/A", "n/a")) %>%
  dplyr::rename_all(snakecase::to_snake_case) %>%
  dplyr::select(-uwbc_batch) %>%
  dplyr::filter(uwbc_submission == "March 2021")


message("separating by batch")

sample_table %<>%
  dplyr::mutate(
    filename = str_split(basename(end1), regex("_R[1|2]")),
    filename = map_chr(filename, 1)) %>%
  dplyr::select(-tidyselect::starts_with("end")) %>%
  dplyr::group_split(batch)

## get old data
wisc_dir <- "/z/Comp/onglab/Projects/WISC"
folder <- "2020_03_10_final_dust_processing"
project <- "dust2M"

data_dir <- file.path(wisc_dir, "data", "microbiome", "house_dust")

old_sample_table <- read_csv(file.path(data_dir, "sample_tables",
  "wisc_dust_microbiome_sample_list_corrected.csv"))



meta <- old_sample_table %>%
  rename(batch = uwbc_batch) %>%
  select(name, subject_id, location, batch) %>%
  rename(baby_sid = subject_id) %>%
  mutate(location = snakecase::to_snake_case(location))

## manual cleaning
meta %<>%
  mutate(
    location = if_else(location == "farm_hsdust",
      "hsdust", location),
    location = if_else(str_detect(name, "Blank") |
      str_detect(name, "Neg"), NA_character_, location),
    batch = str_split(name, "\\_"),
    batch = map_chr(batch, 1))


message("correcting meta")
meta_list <- meta %>%
  group_split(batch)

message("processing December 2018")
meta_dec_2018 <- meta_list[[2]] %>%
  dplyr::select(-batch) %>%
  dplyr::mutate(filename = str_remove(name, "batch2_")) %>%
  dplyr::inner_join(sample_table[[1]], by = "filename") %>%
  dplyr::select(-name) %>%
  dplyr::select(batch, key, tidyselect::everything())

message("processing December 2019")
meta_dec_2019 <- meta_list[[3]] %>%
  dplyr::select(-batch) %>%
  dplyr::mutate(filename = str_remove(name, "batch3_")) %>%
  dplyr::right_join(sample_table[[2]], by = "filename") %>%
  dplyr::select(-name) %>%
  dplyr::select(batch, key, tidyselect::everything()) %>%
  dplyr::filter(!str_detect(filename, regex("^PDcomp"))) %>%
  dplyr::filter(!str_detect(filename, regex("^CTAB")))

message("processing June 2021")
meta_jun_2021 <- xls_table %>%
  dplyr::select(filename, subject_id, location) %>%
  dplyr::right_join(sample_table[[3]], by = "filename") %>%
  dplyr::rename(baby_sid = subject_id) %>%
  dplyr::mutate(location = snakecase::to_snake_case(location)) %>%
  dplyr::select(batch, key, baby_sid, location, filename)
  
meta <- bind_rows(meta_dec_2018, meta_dec_2019, meta_jun_2021)



## add latent classes
latent_classes <- file.path(wisc_dir, "results", "aklobo",
  "2020_03_25_rerun_wisc_lca/new_farm_exposure_classes.csv") %>%
  read_csv()

latent_classes %<>%
  select(class, Baby_SID) %>%
  mutate(Baby_SID = as.character(Baby_SID))

meta %<>%
  left_join(latent_classes, by = c(baby_sid = "Baby_SID"))

## add survey info
redcap_wisc <- file.path(wisc_dir, "extracted_data",
  "redcap", "txt")
redcap_wfs <- file.path(wisc_dir, "extracted_data",
  "amish_redcap", "txt")

### survey: Dmgrphcs
dmgrphcs <- bind_rows(
  file.path(redcap_wisc, "Dmgrphcs.data.2020_11_03.txt") %>%
    read_tsv(),
  file.path(redcap_wfs, "Dmgrphcs.data.2020_09_17.txt") %>%
    read_tsv())

dmgrphcs %<>%
  rename_all(list(snakecase::to_snake_case)) %>%
  select(-birth_season, -pilot,
    -basic_bio_bank,
    -family_id) %>%
  rename(birth_month_cat = birthmonthcat)

dmgrphcs %<>%
  mutate_if(is.character, list(snakecase::to_snake_case)) %>%
  mutate_at(vars(ends_with("sid")), list(as.character))

### survey: ICRA
icra <- bind_rows(
  file.path(redcap_wisc, "ICRA.data.basic.2020_03_13.txt") %>%
    read_tsv(),
  file.path(redcap_wfs, "ICRA.data.basic.2020_09_17.txt") %>%
    read_tsv()
)

icra %<>%
  rename_all(list(snakecase::to_snake_case)) %>%
  mutate_if(is.character, list(snakecase::to_snake_case)) %>%
  mutate_at(vars(ends_with("sid")), list(as.character)) %>%
  select(-icra_wfs_vaccinations)


### survey: PRE

pre <- bind_rows(
  file.path(redcap_wisc, "PRE.data.plus.2020_09_17.txt") %>%
    read_tsv(),
  file.path(redcap_wfs, "PRE.data.plus.2020_09_17.txt") %>%
    read_tsv()
)

pre %<>%
  rename_all(list(snakecase::to_snake_case)) %>%
  mutate_if(is.character, list(snakecase::to_snake_case)) %>%
  mutate_at(vars(ends_with("sid")), list(as.character)) %>%
  select(-months)

meta %<>%
  dplyr::mutate(
    baby_sid = if_else(str_detect(baby_sid, regex("0$")),
      as.character(as.numeric(baby_sid) + 1), baby_sid))


out <- tribble(
  ~ prefix, ~ content, ~ data,
  "meta", "meta information of filenames", meta,
  "dmgrphcs", "demographic information", dmgrphcs,
  "icra", "information related to respiratory data", icra,
  "pre", "pre-birth survey data", pre)

fs::dir_create(dirname(arguments$outfile))
qs::qsave(out, arguments$outfile)
