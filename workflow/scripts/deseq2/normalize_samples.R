source("renv/activate.R")

info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

library(magrittr)
library(phyloseq)
library(DESeq2)
library(qs)

### define directories
message("loading data")

ps <- qs::qread(snakemake@input[["phyloseq"]])


deseq <- phyloseq::phyloseq_to_deseq2(ps, design = ~ 1)
message("removing empty samples (if there is any)")
deseq <- deseq[, colSums(assay(deseq)) > 0] # removing empty samples if any
deseq <- DESeq2::estimateSizeFactors(deseq, type = "poscounts")
deseq_norm <- BiocGenerics::counts(deseq, normalized = TRUE)

ps_norm <- phyloseq::phyloseq(
  otu_table(deseq_norm, taxa_are_rows = TRUE),
  tax_table(ps),
  sample_data(ps))

qs::qsave(ps_norm, snakemake@output[["phyloseq_norm"]])