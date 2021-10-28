#!/ua/chasman/miniconda3/envs/microbiome/bin/Rscript

#' Generates metadata and negative control files
#' specifically for the 2018 and 2019 batches of stool and vaginal swab 16S
#' @author rwelch
#' @author chasman

"Prepares metadata table and negative control mapping table.

Usage:
prepare_metadata_stool.R [options] [--mapping <mapping_table>... --table <input>]
prepare_metadata_stool.R (-h|--help)
prepare_metadata_stool.R --version

Options:
--out_meta=<out_metadata>    Name for output metadata filename [default: data/meta.qs]
--out_control=<out_control>  Name for output negative control mapping file [default: data/negcontrol.qs]
--out_reps=<out_reps>        Name for output replicate table [default: data/replicates.tsv]
--samples=<sample_table>     Name of the file with the samples table [default: samples.tsv]
--mapping=<mapping_table>    Name of file with sample to plate mapping info.
--table=<input_metadata>     Name of an input metadata table file to merge.
--in_data=<data>             Name of input metadata file with R1/R2 filenames [default: data/sample_info.csv]" -> 
  doc
  
library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

if (interactive()) {
  my_args=c(my_args, c("--mapping", "data/Nextseq_04202018_mapping_WISC.txt",
    "--mapping", "data/Nextseq_190826_mapping_WISC.txt",
    "--table", "data/big_metadata.qs"))
}

arguments <- docopt(doc, args = my_args,
  version = "prepare metadata v1")

stopifnot(file.exists(arguments$samples))

library(magrittr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(readxl, quietly = TRUE)
library(qs, quietly = TRUE)

#' reads a tabular file or qs or rds
read_file <- function(filename) {
  ext <- tools::file_ext(filename)
  if (ext %in% c("txt", "tsv")) {
    readr::read_tsv(filename) 
  } else if (ext == "qs") {
    qs::qread(filename)
  } else if (ext == "csv") {
    readr::read_csv(filename)
  } else if (ext == "rds") {
    readRDS(filename)
  } else {
    message("no extension matched")
    NULL
  }
}

sample_table <- read_file(arguments$samples)
sample_table %<>%
  dplyr::mutate(across(c(end1, end2), basename, .names = "{.col}_base")) %>%
  tidyr::separate(end1_base, into = c("sample_name"), sep = "_R", remove = F)

get_sampletype <- function(seq_id) {
  if (grepl("^V|VGSWAB", seq_id)) {
    "VGSWAB00"
  } else if (grepl("STOOL12", seq_id)) {
    "STOOL09_12"
  } else if (grepl("STOOL02", seq_id)) {
    "STOOL02"
  } else if (grepl("MECON", seq_id)) {
    "MECON00"
  } else if (grepl("NEC|AMIES|NTC", seq_id)) {
    "NEC"
  } else if (grepl("^S", seq_id) && !grepl("STOOL02",seq_id)) {
    "STOOL09_12"
  } else {
    "HELP"
  }
}

mapping_tb <- purrr::map(arguments$mapping, 
  ~ read_file(.x) %>% dplyr::mutate(mapping = .x)) %>%  
  bind_rows() %>%
  dplyr::rename(map_id = `#SampleID`) %>%
  dplyr::mutate(sample_type = purrr::map(
      map_id, get_sampletype) %>% unlist) 

# keep only things we need for this analysis
mapping_tb %<>%
  dplyr::select(map_id, sample_type, mapping, PrimerPlate, PrimerWell)

# match sample table filenames to sample info and mapping
sample_info <- read_file(arguments$in_data)
sample_info %<>%
  dplyr::mutate(across(c(R1, R2), basename, .names = "{.col}_base")) %>%
  dplyr::transmute(
     subj_id = Subject,
     flag_reads,
     R1_base, 
     seq_id = SeqID) %>%
  distinct()

# start with sample table, then map in mapping table
sample_meta <- sample_table %>%
  dplyr::select(batch, key, end1_base, end2_base) %>%
  tidyr::separate(end1_base, into = c("map_id"), sep="_", remove=F) %>%
  dplyr::mutate(map_id = purrr::map2(map_id, batch,
    ~{if (grepl("STOOL09", .x) && .y=="batch2019") {
        gsub(".STOOL09.01", "", .x)
     } else {
        .x
     }}) %>% unlist) %>%
  left_join(mapping_tb) %>%
  left_join(sample_info, 
    by = c("end1_base" = "R1_base")) 

stopifnot(nrow(sample_meta) == nrow(sample_table))

# get subject IDs, which can be baby or mom sids
sample_meta %<>%
  dplyr::mutate(subj_id = 
    purrr::map2(subj_id, map_id,
      ~{if (is.na(.x)) {
          str_match(.y, "(\\d{4}).*")[,2]
      } else {
        .x
      }} ) %>% unlist)
stopifnot(nrow(sample_meta) == nrow(sample_table))

# At this point we will add in more SID metadata
# there are some mom sids with multiple baby_sids
# so we'll have to make it a list
metas = tibble()
if (!is.null(arguments$table)) {
  metas <- read_file(arguments$table)
  
  # for mom_sids, remove 
  # any baby-specific features.
  mom_meta <- metas %>%
    dplyr::select(-baby_sid, -sex, -matches("shotgun|02")) %>%
    distinct() %>%
    dplyr::mutate(subj_id = mom_sid)
  re_metas <- metas %>%
    mutate(subj_id = baby_sid) %>%
    bind_rows(mom_meta)

  full_sample_meta <- sample_meta %>%
    dplyr::mutate(mom_sid = round(as.integer(subj_id), -1) %>%
      as.character) %>%
    left_join(re_metas)
  stopifnot(nrow(sample_meta) == nrow(full_sample_meta))

  sample_meta <- full_sample_meta
}

# now map negative controls to batch and plate
negcontrol_tb = dplyr::filter(sample_meta, sample_type=="NEC") %>%
  dplyr::select(batch, map_id, key, PrimerPlate) %>%
  dplyr::group_by(batch, PrimerPlate) %>%
  dplyr::summarize(negcontrols = list(map_id), kits = list(key))

# For some reason, batch2018 plate 3 has no negative controls
# this has only a few samples and they are all repeats.
# I'll assign the VGSWAB to negative controls for plate 2
# and the STOOL02 to negative controls for plate 4
fix_plates = sample_meta %>%
  dplyr::filter(PrimerPlate=="PLATE3" & batch=="batch2018") %>%
  dplyr::mutate(PrimerPlate = ifelse(
      sample_type == "VGSWAB00", "PLATE2", "PLATE4")) 

sample_meta %<>%
  dplyr::filter(!(map_id %in% fix_plates$map_id)) %>%
  bind_rows(fix_plates)

# finalize and save negative control map
negcontrol_map = sample_meta %>%
  dplyr::filter(sample_type != "NEC") %>%
  left_join(negcontrol_tb) %>%
  dplyr::select(batch, key, seq_id, negcontrols, kits)
qs::qsave(negcontrol_map, arguments$out_control)


# define replicates
# take the first one as the main sample
# we will merge the other one into it
main_key <- sample_meta %>%
  dplyr::select(key, subj_id, sample_type, seq_id) %>%
  dplyr::group_by(subj_id, sample_type) %>%
  dplyr::slice_min(order_by = seq_id) %>%
  dplyr::transmute(main_key = key, subj_id, sample_type)
rep_tb <- sample_meta %>%
  dplyr::select(key, subj_id, sample_type) %>%
  left_join(main_key)
qs::qsave(rep_tb, arguments$out_reps)


#################
# save metadata table
sample_meta %>%
  #dplyr::select(-seq_id, -mapping) #%>%
  qs::qsave(arguments$out_meta)

# output tsv
if (grepl("\\.qs$", arguments$out_meta)) {
  out_txt = gsub("\\.qs", ".tsv", arguments$out_meta)
  sample_meta %>%
    #dplyr::select(-seq_id, -mapping) %>%
    write_tsv(out_txt)
}
