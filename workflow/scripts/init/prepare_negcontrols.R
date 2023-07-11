
library(magrittr)
library(tidyverse)

samples <- readr::read_tsv(snakemake@input[["samples"]])
meta <- qs::qread(snakemake@input[["meta"]])

# split samples by batch
samples %<>%
  inner_join(meta, by = c("batch", "key")) %>%
  group_split(batch) %>%
  set_names(c("2018", "2020", "2022"))

# init neg_controls output
neg_controls <- list()

# we are using this function to match the negative controls
# in terms of either unique plate or unique well
match_negative_controls <- function(plate, well, negs) {

  negs1 <- negs %>%
    filter(primer_plate == plate) %>%
    pull(key)
  negs2 <- negs %>%
    filter(primer_well == well) %>%
    pull(key)

  unique(c(negs1, negs2))

}


# 2018
sample_map <- snakemake@params[["batch_2018"]] %>%
  readr::read_tsv() %>%
  rename_with(snakecase::to_snake_case)

sample_map <- samples[["2018"]] %>%
  mutate(
    sample_id = end1 %>%
      basename() %>%
      str_remove("_R1.fastq.gz")) %>%
  left_join(sample_map, by = "sample_id") %>%
  select(-end1, -end2, -contains("sequence"),
    -starts_with("sample"), -owner, -addto_pool, -description)

negs <- sample_map %>%
  filter(is_neg) %>%
  select(key, starts_with("primer"))

neg_controls[["2018"]] <- sample_map %>%
  filter(!is_neg) %>%
  mutate(
    kits = map2(primer_plate, primer_well,
      match_negative_controls, negs)) %>%
  select(batch, key, kits)

# 2020
sample_map <- snakemake@params[["batch_2020"]] %>%
  readr::read_tsv() %>%
  rename_with(snakecase::to_snake_case)

sample_map1 <- samples[["2020"]] %>%
  mutate(
    sample_id = end1 %>%
      basename() %>%
      str_remove("_R1.fastq.gz")) %>%
  inner_join(sample_map, by = "sample_id") %>%
  select(-end1, -end2, -contains("sequence"),
    -starts_with("sample"), -owner, -description)

sample_map2 <- samples[["2020"]] %>%
  mutate(
    sample_id = end1 %>%
      basename() %>%
      str_remove(".STOOL09.01_R1.fastq.gz")) %>%
  inner_join(sample_map, by = "sample_id") %>%
  select(-end1, -end2, -contains("sequence"),
    -starts_with("sample"), -owner, -description)

sample_map <- bind_rows(sample_map1, sample_map2)

negs <- sample_map %>%
  filter(is_neg) %>%
  select(key, starts_with("primer"))

neg_controls[["2020"]] <- sample_map %>%
  filter(!is_neg) %>%
  mutate(
    kits = map2(primer_plate, primer_well,
      match_negative_controls, negs)) %>%
  select(batch, key, kits)

# 2022
sample_map <- snakemake@params[["batch_2022"]] %>%
  readr::read_tsv() %>%
  rename_with(snakecase::to_snake_case) %>%
  filter(sample_id != "S2.4031")

sample_map <- samples[["2022"]] %>%
  mutate(
    sample_id = end1 %>%
      basename() %>%
      str_remove("_R1.fastq.gz")) %>%
  left_join(sample_map, by = "sample_id") %>%
  select(-end1, -end2, -contains("sequence"),
    -starts_with("sample"), -description)

negs <- sample_map %>%
  filter(is_neg) %>%
  select(key, starts_with("primer"))

neg_controls[["2022"]] <- sample_map %>%
  filter(!is_neg) %>%
  mutate(
    kits = map2(primer_plate, primer_well,
      match_negative_controls, negs)) %>%
  select(batch, key, kits)

neg_controls %<>%
  bind_rows()

neg_controls %>%
  qs::qsave(snakemake@output[[1]])
