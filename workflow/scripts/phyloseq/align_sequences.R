#library(dada2)
source("renv/activate.R")

sink(snakemake@log[[1]])
info <- Sys.info();

message(stringr::str_c(names(info), " : ", info, "\n"))
print(stringr::str_c(names(info), " : ", info, "\n"))

library(magrittr)
library(DECIPHER)
library(phangorn)
library(qs)
library(ape)

### define directories
message("loading data")
asv_table <- qs::qread(snakemake@input[["asv"]])

seqs <- colnames(asv_table)
names(seqs) <- stringr::str_c("asv", seq_along(seqs), sep = "_")

#seqs <- getSequences(seqtab)
message("aligning sequences")
alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA)
Biostrings::writeXStringSet(alignment, snakemake@output[["alignment"]])


phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- phangorn::dist.ml(phang.align)
tree_nj <- phangorn::NJ(dm) # Note, tip order != sequence order
tree_nj <- phangorn::midpoint(tree_nj)


ape::write.tree(tree_nj, snakemake@output[["tree"]])


# fit <- phangorn::pml(treeNJ, data = phang.align)

# ## negative edges length changed to 0!


# fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")


# fitGTR <- update(fit, k = 4, inv = 0.2)
# fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#             rearrangement = "stochastic",
#             control = phangorn::pml.control(trace = 0))

