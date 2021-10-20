#!/usr/bin/Rscript

## DADA2 Implementation 
## Author: Katie McCauley, based on the tutorial by Ben Callahan
## Most-recent update: May 13, 2019

# 05/13/19: Add in code at the beginning to install first install pacman (if missing) and then to install other difficult packages without user input.
# 02/14/20: Modifying to run fully on each run instead of merging several runs. Should, theoretically, output run-specific OTU tables that can then be merged as needed.


rm(list=ls())
args <- commandArgs(TRUE)
# This now checks for basic/difficult packages that need to be installed and installs them if not already done.
if(!"pacman" %in% installed.packages()) install.packages("pacman", repos="https://cran.revolutionanalytics.com/")
if(!"devtools" %in% installed.packages()) {
install.packages("devtools")
}
if(!"dada2" %in% installed.packages()) {
devtools::install_github("benjjneb/dada2", ref="v1.12") # change the ref argument to get other versions
}
if(!"DECIPHER" %in% installed.packages()) {
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER")
}

.libPaths()

pacman::p_load(dada2, parallel, purrr, tictoc, DECIPHER, dplyr, tibble)

.libPaths()

source("DADA2_modfunc.R") # I am modifying some of the DADA2 functions. To see specific details, view the "DADA2_modfunc.R" script and search for "(KM)". I have noted my modifications and my reasoning for the change(s). Keep this in the seq_path directory.
tic("Completed DADA2 Pipeline") # This function gives you the total runtime of the pipeline until the "toc" function at the very end. New "tic" functions are started within this that have "toc" functions also within. 
threads <- as.numeric(args[1]) #Number of cores you want to use throughout the process. Modified here to use the number of cores specified in the bash script.
seq_path <- "." # This says that your sequences (and therefore your CSV file and NextSeq directories) are contained within the working directory.
#You can also designate a "save path" (so that you can save your final files in a separate place from your sequence files)
save_path <- "dada2_output_files" # A unique name for saving some auxilary files as well as your output. Can be changed if you want to, say, run with a different subset (based on pulling out samples from your CSV file.
setwd(seq_path)
threshold <- 50 # Drop a sample if it contains less than this threshold of reads (shouldn't be lower than 20).


#This sequence file can either be made by hand or made using my script titled "Code_for_SampListFile.R", and has three columns: one named "fastqs" that has the location of all sample fastqs, one named "sample.name" with the sample names for each of your sample files, and one named "direction" indicating "forward" or "reverse" (ie, the direction of the read file").
csv_file <- "fastq_file_list" # This file name will be used as the downstream output names as well (so that if you, say, run this script on a different cSV file but in the same location, you won't accidentally over-write your work -- not that it's ever happened before...)


###################################################
#####       START DADA2 CODE            ###########
###################################################

if(!is.null(csv_file)) {
file_list <- read.csv(paste0(csv_file, ".csv"),stringsAsFactors=FALSE)
fnFs <- list()
fnRs <- list()
filtFs <- list()
filtRs <- list()
for(i in 1:length(unique(file_list$run.source))) {
fnFs[[i]] <- file_list$fastqs[file_list$direction %in% "R1" & file_list$run.source %in% unique(file_list$run.source)[i]]
fnRs[[i]] <- file_list$fastqs[file_list$direction %in% "R2" & file_list$run.source %in% unique(file_list$run.source)[i]]

filt_path <- file.path(save_path, "filt_fastqs") # Place filtered fastqs in filtered/ subdirectory

filtFs[[i]] <- file.path(filt_path, paste0(file_list$sample.names[file_list$direction %in% "R1" & file_list$run.source %in% unique(file_list$run.source)[i]], "_R1_filt.fastq.gz"))
filtRs[[i]] <- file.path(filt_path, paste0(file_list$sample.names[file_list$direction %in% "R2" & file_list$run.source %in% unique(file_list$run.source)[i]], "_R2_filt.fastq.gz"))
}

} else {
# This part isn't well developed anymore, but it is used if there isn't a .csv file with the run information. Maybe make this part obsolete? Willing to further develop if desired.
fnFs <- sort(list.files(seq_path, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(seq_path, pattern="R2_001.fastq.gz", full.names = TRUE))

#Ensure that your forward and reverse read file paths are correctly 
print(fnFs)
print(fnRs)

#I'll use the defaults that we currently employ, given that it's too much work to do this sample-by-sample?
filt_path <- file.path(path, "filt_fastqs") # Place filtered files in filtered/ subdirectory
sample.names <- sapply(strsplit(basename(fnFs), "[_]"), function(x) x[1])# My sample names were kinda complicated, so I needed to extract the sample names and create cleaner sample name files for those that get filtered.
filtFs <- file.path(filt_path, paste0(sample.names, "_fwd_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_rev_filt.fastq.gz"))
}


for(i in 1:length(fnFs)) {
print("Begin the Filter and Trim Step")
out <- filterAndTrim2(fwd=unlist(fnFs[[i]]), filt=unlist(filtFs[[i]]), rev=unlist(fnRs[[i]]), filt.rev=unlist(filtRs[[i]]), 
              maxEE=2, truncQ=0, rm.phix=TRUE, minLen=150,
              compress=TRUE, multithread=threads, verbose=TRUE, matchIDs = TRUE)

prop.out <- out[,"reads.out"]/out[,"reads.in"]
out <- cbind(out, prop.out)
print(out)

keep <- out[,"reads.out"] > threshold

drop.names <- gsub("_R1.fastq.gz", "", names(keep[!keep]))
if(length(drop.names)>0) {

filtFs[[i]] <- filtFs[[i]][!grepl(paste(drop.names, collapse="|"), filtFs[[i]])]
filtRs[[i]] <- filtRs[[i]][!grepl(paste(drop.names, collapse="|"), filtRs[[i]])]

print(paste("Samples were dropped due to having fewer than", threshold, "reads per sample:"))
print(paste(drop.names, collapse=","))
}

tic("Complete Procesing Up To mergePairs()")
derepFs <- list()
derepRs <- list()
errF <- list()
errR <- list()
dadaFs <- list()
dadaRs <- list()

## Dereplication
print("Dereplicate Sequences")
tic("Dereplication")
derepFs[[i]] <- derepFastq(filtFs[[i]], verbose=TRUE)
derepRs[[i]] <- derepFastq(filtRs[[i]], verbose=TRUE)
toc()
tic("Error-learning")
errF[[i]] <- learnErrors(derepFs[[i]], multithread=threads, randomize=TRUE, nbases=1e8)
errR[[i]] <- learnErrors(derepRs[[i]], multithread=threads, randomize=TRUE, nbases=1e8)
toc()
## Run DADA2
tic("DADA2 Step")
print("DADA2 Step")
dadaFs[[i]] <- dada(derepFs[[i]], err=errF[[i]], pool=FALSE, multithread=threads)
dadaRs[[i]] <- dada(derepRs[[i]], err=errR[[i]], pool=FALSE, multithread=threads)
toc()
 # This will be the end of the run-specific processing for now.
toc()
#save.image(paste0(save_path, "/dada2_result_object_postdada.RData"))

## NOW Merge Paired Ends
print("Merge Paired Reads")
library(purrr) #bringing in this package because it "undoes" the list-of-lists that I had, which are no longer needed because I have already taken care of the run-specific processing
mergers <- mergePairs(flatten(dadaFs), flatten(derepFs), flatten(dadaRs), flatten(derepRs), minOverlap=25, verbose=TRUE)

save.image(paste0(save_path, "/dada2_result_object.RData"))

#I think this is useful information to print, so I will keep it here.
head(mergers[[1]])
str(mergers[[1]])
dim(mergers[[2]])

#Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
sum(seqtab)

#A data check..
print("Initial Sequence Lengths Pre-Chimera Removal")
table(nchar(getSequences(seqtab)))
# Distribution of sequence reads:
seq.length <- sort(table(nchar(getSequences(seqtab))))
## Remove Chimeras:

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
        multithread=threads, verbose=TRUE)
dim(seqtab.nochim)
print("Proportion of non-chimeric reads")
sum(seqtab.nochim)/sum(seqtab)

print("Initial Sequence Lengths")
table(nchar(getSequences(seqtab.nochim))) # Moving this to *after* the bimera removal, since it's possible that the odd-sized reads are considered bimeras.
seq.length <- as.numeric(names(which.max(table(nchar(getSequences(seqtab.nochim))))))
seq.ideal <- seq.int(from=seq.length-5, to=seq.length+5)
seqtab.nochim <- seqtab.nochim[, nchar(getSequences(seqtab.nochim)) %in% seq.ideal]

print("New Sequence Length Distribution")
table(nchar(getSequences(seqtab.nochim)))

save.image(paste0(save_path, "/dada2_result_object.RData"))

getN <- function(x) sum(getUniques(x))
sumdadaFs <- as.data.frame(sapply(flatten(dadaFs), getN)) %>%
	rownames_to_column(var="sampleid") %>%
	mutate(sampleid = gsub( "_R1_filt.fastq.gz","",sampleid)) %>%
	rename('denoisedF' = 'sapply(flatten(dadaFs), getN)')
sumdadaRs <- as.data.frame(sapply(flatten(dadaRs), getN)) %>%
        rownames_to_column(var="sampleid") %>%
        mutate(sampleid = gsub("_R2_filt.fastq.gz","",sampleid)) %>%
        rename('denoisedR' = 'sapply(flatten(dadaRs), getN)')
sum.merge <- as.data.frame(sapply(mergers, getN)) %>%
        rownames_to_column(var="sampleid") %>%
        mutate(sampleid = gsub("_R1_filt.fastq.gz","",sampleid)) %>%
        rename('merged' = 'sapply(mergers, getN)')
sum.seqtab <- as.data.frame(rowSums(seqtab.nochim)) %>%
        rownames_to_column(var="sampleid") %>%
        mutate(sampleid = gsub("_R1_filt.fastq.gz","",sampleid)) %>%
        rename('nonchim' = 'rowSums(seqtab.nochim)')
sum.out <- as.data.frame(out) %>%
	rownames_to_column(var="sampleid") %>%
	mutate(sampleid = gsub("_R1.fastq.gz","", sampleid))
track <- sum.out %>%
	full_join(sumdadaFs) %>%
	full_join(sumdadaRs) %>%
	full_join(sum.merge) %>%
	full_join(sum.seqtab)

write.csv(track, paste0(save_path, "/TrackedReadsThruDADA2.csv"))
rm(list=c("dadaFs","dadaRs","mergers"))

## Taxonomy

#Need to figure out where Silva will end up going. Right now, I have copied it into the directory for this run, but that's not sustainable.
#Discussion of whether this should be DECIPHER or the usual taxonomy assignment. I will generate both and leave it up to the user.
print("Running the assignTaxonomy taxonomy Assignment")
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=threads, tryRC=TRUE, minBoot=80,outputBootstraps=TRUE)
write.csv(taxa$boot, paste0(save_path, "/TaxonomyBootstraps.csv"), row.names=F)
taxa <- addSpecies(taxa$tax, "silva_species_assignment_v132.fa.gz", allowMultiple=TRUE, verbose=TRUE,tryRC=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
table(is.na(taxa.print[,"Genus"]))

print("Running DECIPHER taxonomy assignment method")
dna <- DNAStringSet(getSequences(seqtab.nochim))
load("SILVA_SSU_r132_March2018.RData") # this takes awhile
ids <- IdTaxa(dna, trainingSet, strand="both", processors=threads, verbose=TRUE, threshold=60)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        conf <- x$confidence[max(m, na.rm=TRUE)]
        print(conf)
        taxa.d <- x$taxon[m] #taxa.d is used because I don't want to over-write the 'taxa' object from above.
        taxa.d[startsWith(taxa.d, "unclassified_")] <- NA
        c(taxa.d, conf)
}))
colnames(taxid) <- c(ranks, "confidence")


# Add back "nice" sample names 
if(length(sum.out$sampleid) %in% length(sumdadaRs$sampleid))  rownames(seqtab.nochim) <- sumdadaRs$sampleid

write.table(x=getSequences(seqtab.nochim), file=paste0(save_path,"/", csv_file,"_",unique(file_list$run.date)[i],"_seqs.txt"), sep="\t", col.names=F, row.names=F, quote=F)
write.table(x=t(seqtab.nochim), file=paste0(save_path,"/",csv_file,"_",unique(file_list$run.date)[i],"_otutable.txt"), sep="\t", col.names=T, row.names=T, quote=F)
write.table(x=taxa.print, file=paste0(save_path,"/",csv_file,"_",unique(file_list$run.date)[i], "_tax_table.txt"), sep="\t", col.names=T, row.names=T, quote=F)
write.table(x=taxid, file=paste0(save_path, "/", csv_file,"_",unique(file_list$run.date)[i], "_tax_table_DECIPHER.txt"), sep="\t", col.names=T, row.names=T, quote=F)
}
save.image("dada2_result_object.RData")
toc()
