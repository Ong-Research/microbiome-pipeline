source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(phyloseq)
library(doParallel)
library(foreach)
library(parallelDist)
library(qs)

snakemake@config[["beta"]] <- stringr::str_to_lower(snakemake@config[["beta"]])

message("reading phyloseq from ", snakemake@input[["phyloseq"]])

ps <- qs::qread(snakemake@input[["phyloseq"]])

metrics <- c("bray", "fJaccard", "euclidean",
  "hellinger", "mahalanobis", "manhattan", "bhjattacharyya", "canberra",
  "chord", "unifrac", "w_unifrac", "w_unifrac_norm")
stopifnot(snakemake@config[["beta"]] %in% metrics)

message("computing beta diversity with ", snakemake@config[["beta"]])

if (snakemake@config[["beta"]] %in% c("bray", "fJaccard", "euclidean",
  "hellinger", "mahalanobis", "manhattan", "bhjattacharyya", "canberra",
  "chord")) {

    cl <- parallel::makeCluster(snakemake@threads)
    doParallel::registerDoParallel(cl)

    asv_table <- phyloseq::otu_table(ps, taxa_are_rows = FALSE)
    asv_table <- as.matrix(asv_table@.Data)
    distance <- parallelDist::parDist(asv_table,
      method = snakemake@config[["beta"]], threads = snakemake@threads)

    pcoa <- ape::pcoa(distance)
# distances <- list()

# distances[["unifrac"]] <- UniFrac(asv_ps, weighted = FALSE,
# 	normalized = FALSE, parallel = TRUE)
# distances[["w_unifrac"]] <- UniFrac(asv_ps, weighted = TRUE,
# 	normalized = FALSE, parallel = TRUE)
# distances[["w_unifrac_norm"]] <- UniFrac(asv_ps, weighted = TRUE,
# 	normalized = TRUE, parallel = TRUE)

} else if (snakemake@config[["beta"]] == "unifrac") {

  distance <- phyloseq::UniFrac(ps, weighted = FALSE, normalized = FALSE,
    parallel = TRUE)
  pcoa <- ape::pcoa(distance)

} else if (snakemake@config[["beta"]] == "w_unifrac") {

  distance <- phyloseq::UniFrac(ps, weighted = TRUE, normalized = FALSE,
    parallel = TRUE)
  pcoa <- ape::pcoa(distance)

} else if (snakemake@config[["beta"]] == "w_unifrac_norm") {

  distance <- phyloseq::UniFrac(ps, weighted = TRUE, normalized = TRUE,
    parallel = TRUE)
  pcoa <- ape::pcoa(distance)


}

distance %>% qs::qsave(snakemake@output[["distance"]])
pcoa %>% qs::qsave(snakemake@output[["pcoa"]])
