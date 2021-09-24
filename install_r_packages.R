options("repos" =
  c(CRAN = "https://repo.miserver.it.umich.edu/cran/"))

# need to install libgit2 before

# vscode package
install.packages("lintr")
install.packages("languageserver")
install.packages("usethis")
install.packages("devtools")
install.packages("here")
install.packages("snakecase")

# tidyverse and use
install.packages("magrittr")
install.packages("tibble")
install.packages("tidyverse")
install.packages("tidygraph")
install.packages("tidymodels")
install.packages("stacks")
install.packages("pdp")
install.packages("vip")
install.packages("spelling")
install.packages("sessioninfo")

# plotting
install.packages("ggplot2")
install.packages("scales")
install.packages("ggrepel")
install.packages("pals")
install.packages("cowplot")
install.packages("patchwork")
install.packages("ggraph")

# tables
install.packages("broom")
install.packages("gt")
install.packages("gtsummary")

# read files
install.packages("qs")
install.packages("jsonlite")
install.packages("yaml")

# Rcpp
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("RcppEigen")

# scripting
install.packages("docopt")

# extra
install.packages("taxizedb")

# documentation
install.packages("gert")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("workflowr")
install.packages("Cairo")
install.packages("magick")

# Bioconductor stuff
install.packages("BiocManager")
BiocManager::install("S4Vectors")
BiocManager::install("IRanges")
BiocManager::install("Biostrings")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicAlignments")
BiocManager::install("BiocParallel")
BiocManager::install("dada2")
BiocManager::install("ggtree")
BiocManager::install("phyloseq")
BiocManager::install("DirichletMultinomial")
BiocManager::install("microbiome")
BiocManager::install("mia")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("DECIPHER")
BiocManager::install("decontam")
BiocManager::install("mixOmics")
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")
