options("repos" =
  c(CRAN = "https://repo.miserver.it.umich.edu/cran/"))

# need to install libgit2 before
install.packages("BiocManager")


# vscode package
BiocManager::install("lintr", update = FALSE)
BiocManager::install("languageserver", update = FALSE)
BiocManager::install("usethis", update = FALSE)
BiocManager::install("devtools", update = FALSE)
BiocManager::install("here", update = FALSE)
BiocManager::install("snakecase", update = FALSE)

# tidyverse and use
BiocManager::install("magrittr", update = FALSE)
BiocManager::install("tibble", update = FALSE)
BiocManager::install("tidyverse", update = FALSE)
BiocManager::install("tidygraph", update = FALSE)
BiocManager::install("tidymodels", update = FALSE)
BiocManager::install("stacks", update = FALSE)
BiocManager::install("pdp", update = FALSE)
BiocManager::install("vip", update = FALSE)
BiocManager::install("spelling", update = FALSE)
BiocManager::install("sessioninfo", update = FALSE)

# plotting
BiocManager::install("ggplot2", update = FALSE)
BiocManager::install("scales", update = FALSE)
BiocManager::install("ggrepel", update = FALSE)
BiocManager::install("pals", update = FALSE)
BiocManager::install("cowplot", update = FALSE)
BiocManager::install("patchwork", update = FALSE)
BiocManager::install("ggraph", update = FALSE)
BiocManager::install("ggdist", update = FALSE)
BiocManager::install("ggtext", update = FALSE)
BiocManager::install("tidytext", update = FALSE)
devtools::install_github("teunbrand/ggh4x")

# tables
BiocManager::install("broom", update = FALSE)
BiocManager::install("gt", update = FALSE)
BiocManager::install("gtsummary", update = FALSE)

# read files
BiocManager::install("qs", update = FALSE)
BiocManager::install("jsonlite", update = FALSE)
BiocManager::install("yaml", update = FALSE)

# Rcpp
BiocManager::install("Rcpp", update = FALSE)
BiocManager::install("RcppArmadillo", update = FALSE)
BiocManager::install("RcppEigen", update = FALSE)

# scripting
BiocManager::install("docopt", update = FALSE)

# extra
BiocManager::install("parallelDist", update = FALSE)
BiocManager::install("taxizedb", update = FALSE)
BiocManager::install("zCompositions", update = FALSE)
BiocManager::install("Hmisc", update = FALSE)

# documentation
BiocManager::install("gert", update = FALSE)
BiocManager::install("knitr", update = FALSE)
BiocManager::install("rmarkdown", update = FALSE)
BiocManager::install("workflowr", update = FALSE)
BiocManager::install("Cairo", update = FALSE)
BiocManager::install("magick", update = FALSE)

# Bioconductor stuff
BiocManager::install("S4Vectors", update = FALSE)
BiocManager::install("IRanges", update = FALSE)
BiocManager::install("Biostrings", update = FALSE)
BiocManager::install("GenomicRanges", update = FALSE)
BiocManager::install("GenomicAlignments", update = FALSE)
BiocManager::install("BiocParallel", update = FALSE)
BiocManager::install("dada2", update = FALSE)
BiocManager::install("ggtree", update = FALSE)
BiocManager::install("phyloseq", update = FALSE)
BiocManager::install("corncob", update = FALSE)
BiocManager::install("DirichletMultinomial", update = FALSE)
BiocManager::install("microbiome", update = FALSE)
BiocManager::install("mia", update = FALSE)
BiocManager::install("SummarizedExperiment", update = FALSE)
BiocManager::install("scater", update = FALSE)
BiocManager::install("DESeq2", update = FALSE)
BiocManager::install("DECIPHER", update = FALSE)
BiocManager::install("decontam", update = FALSE)
BiocManager::install("mixOmics", update = FALSE)
BiocManager::install("circlize", update = FALSE)
BiocManager::install("ComplexHeatmap", update = FALSE)
BiocManager::install("ComplexUpset", update = FALSE)
BiocManager::install("stringdist", update = FALSE)
BiocManager::install("fastreeR", update = FALSE)
