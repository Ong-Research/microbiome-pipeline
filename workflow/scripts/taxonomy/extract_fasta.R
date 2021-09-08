source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(magrittr)
library(tidyverse)
library(Biostrings)
library(qs)

asvs <- qs::qread(snakemake@input[["asv"]])
sequences <- colnames(asvs)
names(sequences) <- stringr::str_c("asv", seq_along(sequences), sep = "_")


sequences <- Biostrings::DNAStringSet(sequences)
Biostrings::writeXStringSet(sequences, filepath = snakemake@output[["fasta"]])
