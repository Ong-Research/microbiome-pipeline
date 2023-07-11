
# example of how the metadata was cleaned for stool samples
# ---- edit below ----

library(magrittr)
library(tidyverse)

samples <- read_tsv(snakemake@input[[1]])
samples %<>% group_split(batch) %>%
  set_names(c("2018", "2020", "2022"))


# manual curation
meta <- list()

## split by batch

# 2018 batch
meta[["2018"]] <- samples[["2018"]] %>%
  mutate(
    aux = basename(end1),
    is_neg = str_detect(aux, "NTC"),
    months = if_else(str_detect(aux, "STOOL02") & ! is_neg,
      "2M", NA_character_),
    baby_sid = str_split(aux, "\\.") %>% map_chr(1),
    baby_sid = if_else(is_neg, NA_character_, baby_sid),
    mom_sid = as.character(as.numeric(baby_sid) - 1)) %>%
  select(-aux, -end1, -end2) %>%
  select(batch, key, ends_with("sid"), months, is_neg)

# 2020 batch
meta[["2020"]] <- samples[["2020"]] %>%
  mutate(
    aux = basename(end1),
    is_neg = str_detect(aux, "NEC"),
    months = if_else(str_detect(aux, "STOOL02"), "2M", "9M"),
    months = if_else(is_neg, NA_character_, months),
    baby_sid = str_split(aux, "\\.") %>% map_chr(1),
    baby_sid = if_else(is_neg, NA_character_, baby_sid),
    baby_sid = str_remove(baby_sid, "S"),
    mom_sid = as.character(as.numeric(baby_sid) - 1)) %>%
  select(-aux, -end1, -end2) %>%
  select(batch, key, ends_with("sid"), months, is_neg)

# 2022 batch
meta[["2022"]] <- samples[["2022"]] %>%
  mutate(
    aux = basename(end1),
    is_neg = str_detect(aux, "NTC"),
    months = if_else(is_neg, NA_character_, "9M"),
    baby_sid = str_split(aux, "\\.") %>% map_chr(2),
    baby_sid = if_else(is_neg, NA_character_, baby_sid),
    baby_sid = str_remove(baby_sid, "_R1"),
    mom_sid = as.character(as.numeric(baby_sid) - 1)) %>%
  select(-aux, -end1, -end2) %>%
  select(batch, key, ends_with("sid"), months, is_neg)

meta %<>% bind_rows() # 651 samples

meta %>%
  qs::qsave(snakemake@output[[1]])

