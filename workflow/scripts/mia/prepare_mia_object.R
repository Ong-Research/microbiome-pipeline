#!/usr/local/bin/Rscript

#' Prepara mia object
#' @author rwelch2

"Prepare mia object

Usage:
prepare_mia_object.R [<mia_file>] [--asv=<asv_file> --taxa=<taxa_file> --tree=<tree_file> --meta=<meta_file>] [--asv_prefix=<prefix> --log=<logfile> --config=<config> --cores=<cores>]
prepare_mia_object.R (-h|--help)
prepare_mia_object.R --version

Options:
-h --help    show this screen
--asv=<asv_file>    ASV matrix file
--taxa=<taxa_file>    Taxa file
--tree=<tree_file>    Tree file
--meta=<meta_file>    Metadata file
--log=<logfile>    name of the log file [default: logs/filter_and_trim.log]
--config=<config>    name of the config file [default: config/config.yaml]
--cores=<cores>    number of parallel CPUs [default: 8]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "prepare mia file V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {

  arguments$asv <- "output/dada2/remove_chim/asv_mat_wo_chim.qs"
  arguments$taxa <- "output/taxa/kraken_merged/kraken_rdata.qs"
  arguments$tree <- "output/phylotree/newick/tree.nwk"
  arguments$meta <- "output/init/metadata.qs"
  arguments$asv_prefix <- "asv"

}

print(arguments)

info <- Sys.info();

message("loading packages")
library(magrittr)
library(tidyverse)
library(TreeSummarizedExperiment)
library(Biostrings)
library(BiocParallel)
library(ape)
library(mia)
library(qs)
library(scater)
library(yaml)
library(parallelDist)

stopifnot(
  file.exists(arguments$asv),
  file.exists(arguments$taxa),
  file.exists(arguments$tree),
  file.exists(arguments$meta)
)

if (!is.null(arguments$config)) stopifnot(file.exists(arguments$config))

asv <- qs::qread(arguments$asv)
taxa <- qs::qread(arguments$taxa)
tree <- ape::read.tree(arguments$tree)

if (tools::file_ext(arguments$meta) == "tsv") {
  meta <- readr::read_tsv(arguments$meta)
} else if (tools::file_ext(arguments$meta) == "qs") {
  meta <- qs::qread(arguments$meta)
}

cdata <- meta %>%
  as.data.frame() %>%
  tibble::column_to_rownames("key")

# clean samples
sample_names <- intersect(
  rownames(cdata),
  rownames(asv))

# clean ASVs
asv_sequences <- colnames(asv)
colnames(asv) <- str_c(arguments$asv_prefix, seq_along(asv_sequences),
  sep = "_")
names(asv_sequences) <- colnames(asv)
asv_sequences <- Biostrings::DNAStringSet(asv_sequences)

taxa %<>%
  dplyr::add_count(asv) %>%
  dplyr::filter(n == 1) %>%
  dplyr::select(asv, taxa, taxid, domain, phylum, class, order, family,
    genus, species, kraken_db)

asv_names <- intersect(names(asv_sequences), taxa$asv)

asv <- asv[, asv_names]
asv_sequences <- asv_sequences[asv_names, ]


tree$tip.label %<>%
  stringr::str_sub(2, nchar(.) - 1) # either qiime2 or ape is adding "'" at the 
  # start and end of ASV names

tree <- ape::keep.tip(tree, intersect(asv_names, tree$tip.label))

rdata <- taxa %>%
    as.data.frame() %>%
    tibble::column_to_rownames("asv")

out <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(counts = t(asv[sample_names,])),
  colData = cdata[sample_names, ],
  rowData = rdata,
  rowTree = tree)
  
metadata(out)[["date_processed"]] <- Sys.Date()
metadata(out)[["sequences"]] <- asv_sequences

fs::dir_create(dirname(arguments$mia_file))
qs::qsave(out, arguments$mia_file)
message("done!")

# bpp <- BiocParallel::MulticoreParam(workers = as.numeric(arguments$cores))

# message("estimate diversity")
# out <- mia::estimateDiversity(out, abund_values = "counts", BPPARAM = bpp)

# message("estimate richness")
# out <- mia::estimateRichness(out, abund_values = "counts", BPPARAM = bpp)

# if (!is.null(arguments$config)) {

#   config <- yaml::read_yaml(arguments$config)[["beta"]]

# }

# compute_mds_wrap <- function(mia, dist_name, altexp_name, ncomps) {

#   message("computing beta diversity for ", dist_name)

#   if (dist_name %in% c("bray", "euclidean", "hellinger", "mahalanobis",
#     "manhattan", "bhjattacharyya", "canberra", "chord")) {

#     mia <- scater::runMDS(
#       mia, FUN = parallelDist::parDist,
#       method = dist_name, threads = as.numeric(arguments$cores),
#       name = altexp_name, ncomponents = ncomps, exprs_values = "counts",
#       keep_dist = FALSE)

#   } else if (dist_name == "fjaccard") {

#     mia <- scater::runMDS(
#       mia, FUN = parallelDist::parDist,
#       method = "fJaccard", threads = as.numeric(arguments$cores),
#       name = altexp_name, ncomponents = ncomps, exprs_values = "counts",
#       keep_dist = FALSE)

#   } else if (dist_name == "unifrac") {

#     unifrac <- mia::calculateUniFrac(
#       x = t(counts(mia)), tree = rowTree(mia),
#       weighted = FALSE, normalized = FALSE, BPPARAM = bpp)

#     pcoa <- stats::cmdscale(unifrac, k = ncomps)
#     SingleCellExperiment::reducedDim(mia, altexp_name) <- pcoa


#   } else if (dist_name == "w_unifrac") {
#     unifrac <- mia::calculateUniFrac(
#       x = t(counts(mia)), tree = rowTree(mia),
#       weighted = TRUE, normalized = FALSE, BPPARAM = bpp)

#     pcoa <- stats::cmdscale(unifrac, k = ncomps)
#     SingleCellExperiment::reducedDim(mia, altexp_name) <- pcoa
#   } else if (dist_name == "w_unifrac_norm") {

#     unifrac <- mia::calculateUniFrac(
#       x = t(counts(mia)), tree = rowTree(mia),
#       weighted = TRUE, normalized = TRUE, BPPARAM = bpp)
#     pcoa <- stats::cmdscale(unifrac, k = ncomps)
#     SingleCellExperiment::reducedDim(mia, altexp_name) <- pcoa

#   } else {
#     stop(dist_name, " distance not available")
#   }
#   mia

# }



# if (!is.null(config$beta_div)) {
#   for (div in config$beta_div) {
#     out <- compute_mds_wrap(out, div, str_c("all_", div), config$comps)
#   }
# }

# if (any(config$beta_group != "all")) {

#   beta_group <- config$beta_group
#   beta_group <- beta_group[beta_group != "all"]
#   grouped <- purrr::map(beta_group,
#     ~ mia::agglomerateByRank(out, rank = ., na.rm = FALSE))
#   names(grouped) <- beta_group

#   config$beta_div <- config$beta_div[!str_detect(config$beta_div, "unifrac")]

#   for (gg in names(grouped)) {
#   message(gg)
#     for (div in config$beta_div) {
#       grouped[[gg]] <- compute_mds_wrap(grouped[[gg]], div,
#         str_c(gg, div, sep = "_"), config$comps)
#     }
#   }

#   SingleCellExperiment::altExps(out) <- grouped
# }

