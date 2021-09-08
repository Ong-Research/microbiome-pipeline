source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(phyloseq)
library(microbiome)
library(qs)

message("reading phyloseq from ", snakemake@input[["phyloseq"]])
ps <- qs::qread(snakemake@input[["phyloseq"]])


message("computing diversity")
print("computing diversity")


divs <- c("observed", "chao1", "diversity_inverse_simpson",
  "diversity_gini_simpson", "diversity_shannon", "diversity_coverage",
  "rarity_log_modulo_skewness", "rarity_low_abundance", 
  "rarity_rare_abundance", "dominance_dbp", "dominance_dmn",
  "dominance_absolute", "dominance_relative", "dominance_simpson",
  "dominance_core_abundance")


suppressMessages({
diversity_list <- purrr::map(divs,
  ~ microbiome::alpha(ps, .))
})

diversity_list <- purrr::map(diversity_list,
  tibble::as_tibble, rownames = "sample")


diversity <- purrr::reduce(diversity_list,
  purrr::partial(dplyr::inner_join, by = "sample"))

message("saving alpha diversity in ", snakemake@output[["alpha"]])
diversity %>% qs::qsave(snakemake@output[["alpha"]])
