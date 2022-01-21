#!/usr/local/bin/Rscript

#' Parses kraken2 summary file (taxonomic hierarchy) into a table
#' that maps sequence IDs to taxonomy levels.
#' Produces two outputs: the full table (with all hierarchy levels
#' present in the original summary file, which may include things like subspecies)
#' and another table with standardized levels: 
#' domain, phylum, class, order, family, genus, species
#' @param kraken_summary summary report computed by kraken2
#' @author chasman

"Parse kraken2 summary into taxonomy table

Usage:
parse_kraken_summary.R [<outfile_standard> <outfile_full> <kraken_summary>] [--log=<logfile>]
parse_kraken_summary.R (-h|--help)
parse_kraken_summary.R --version

Options:
--log=<logfile>    name of the log file [default: ./parse_kraken_summary.log]" -> doc

library(docopt)

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(doc, args = my_args,
  version = "parse kraken2 summary file V1")

if (!interactive()) {
  log_file <- file(arguments$log, open = "wt")
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

if (interactive()) {
  db = "minikraken"
  arguments$outfile_standard <- "output_test"
  arguments$outfile_full <- "output_full_table"
  arguments$kraken_summary <- sprintf("output/taxa/kraken/%s/kraken_summary.out", db)
}

info <- Sys.info();
print(stringr::str_c(names(info), " : ", info, "\n"))

message("loading packages")
library(tidyverse)
library(magrittr)

stopifnot(file.exists(arguments$kraken_summary))


# for assembling taxonomy string in mpa style
# we'll only take the main headers
# change prefix for domain to k "kingdom"
mpa_levels <- c("D", "P", "C", "O", "F", "G", "S")
new_colname <- c("domain", "phylum", "class", 
    "order", "family", "genus", "species")
mpa_tib <- tibble(colname = mpa_levels, 
            new_colname = new_colname,
            prefix = tolower(mpa_levels)) %>%
        dplyr::mutate(prefix = ifelse(prefix=="d", "k", prefix))


# read summary file
summ_tb <- readr::read_tsv(arguments$kraken_summary,
    col_names = c("frac_asvs", "n_asvs_below", 
        "n_asvs_assigned", "taxlevel", "taxid", 
        "taxname"),
    trim_ws=F)

# Compute number of spaces - used to determine hierarchical relationships
summ_tb %<>%
    dplyr::mutate(nspace = 
        #purrr::map(taxname, ~ str_count(.x, pattern=" "))) %>%
        purrr::map(taxname, ~ str_match(.x, "(^\\s*)")[,1] %>% nchar)) %>%
    unnest(nspace) %>%
    dplyr::mutate(taxid = as.character(taxid)) %>%
    dplyr::mutate(taxname = purrr::map(taxname, str_trim)) %>%
    unnest(taxname)

#' function to process the summary file into a
#' list of taxonomy strings, indexed by taxID.
#' it's a bit sloppy
search_tree <- function(summ_tb, curr=1, prev=0,
    curr_string = list(tmp="tmp"), curr_nspace = c(-1), debug=F) {

    # check to make sure in synch
    stopifnot(length(curr_string) == length(curr_nspace))
    
    cl <- summ_tb[curr,] %>% as.list()

    # previous space level is at the end of the list
    prev_space <- curr_nspace[length(curr_nspace)]

    # the downstream strings that we will collect
    strings = list()

    # stop if we ran off the end
    if (curr > nrow(summ_tb)) {
        #message("ran off the end")
        return(strings)
    }

    if (debug) {
        message("curr ", curr, " prev ", prev)
        message("prev nspace ", prev_space, " curr nspace ", cl$nspace)
        message("\t\tlength check OK?",  length(curr_nspace) == length(curr_string))
    }
    res <- "whoa"
    if (cl$nspace > prev_space) {

        # if we are still proceeding through the tree, keep going.
        curr_string[[cl$taxlevel]] = cl$taxname        
        curr_nspace <- c(curr_nspace, cl$nspace)

        if (debug) {
            message("moving forward")
            message(sprintf("%s_%s", names(curr_string), curr_string) %>%
                paste(collapse="|"))
            message("\t", paste(curr_nspace, collapse = "|"))
        }
    } else if (cl$nspace <= prev_space) { 
        # if the current level is at or above the prev, 
        # then: 
        # back up until we get to the 
        # back up the string until we get to 
        # the level above the current,
        # add current, 
        # set the prev to current
        # and current to current+1
        if (debug) {
            message("backing up. prev was: ", paste(curr_string, collapse = "|"))
            message("\t\tnspace: ", paste(curr_nspace, collapse = "|"))
        }
        #temp_nspace = curr_nspace
        curr_nspace %<>% .[. < cl$nspace ]
        #curr_string %<>% .[temp_nspace < cl$nspace] 
        curr_string %<>% .[1:length(curr_nspace)]

        if (debug) {
            message("now prev is: ", paste(curr_string, collapse = "|"))
            message("\t\tnspace: ", paste(curr_nspace, collapse = "|"))
        }
        curr_string[[cl$taxlevel]] = cl$taxname
        curr_nspace <- c(curr_nspace, cl$nspace)
    }

    if (debug) {
        message(paste(curr_string, collapse = "|"))
        message("\t", paste(curr_nspace, collapse = "|"))
    }

    strings[[cl$taxid]] <- curr_string
    res <- search_tree(summ_tb, 
            curr = curr+1, prev = curr,
            curr_string, curr_nspace, debug = debug)

    return(c(strings, res))
}

# function to pivot tax strings into a tibble
tax_tibble <- function(strings) {
    tmp <- tibble(taxid = names(strings),
        taxlist = strings) %>%
        dplyr::mutate(taxonomy = 
            purrr::map(taxlist, 
                ~ tibble(taxlevel = names(.x), name = .x) %>%
                unnest(name))) %>%
        unnest(taxonomy)
    tmp %>%
        pivot_wider(names_from = taxlevel, values_from = name) 
}

#' function to convert a taxonomy list to an mpa-style string
#' like k__domain|p__phylum|c__class
mpa_string <- function(taxlist, mpa_tib) {
    sub_mpa <- mpa_tib %>%
        dplyr::filter(colname %in% names(taxlist))
    if ("U" %in% names(taxlist) && !is.na(taxlist[["U"]])) {
        taxlist[["U"]]
    } else {
        taxlist[sub_mpa$colname] %>%
            sprintf("%s__%s", sub_mpa$prefix, .) %>%
            paste(collapse="|") %>%
            gsub(" ", "_", .)
    }
}


# Obtain taxonomy strings from table
strings <- summ_tb %>%
    search_tree(debug = T)

# get order of taxlevels
all_levels = c("U", "R", "D", "P", "C", "O", "F", "G", "S")
all_levels = purrr::map(all_levels,
    ~ c(.x, sprintf("%s%s", .x, 1:5))) %>% 
    unlist

# Make table with full taxonomy
# with all levels from original file
# put in order of taxonomy
taxtib <- tax_tibble(strings) %>%
    dplyr::mutate(tax_string = 
        purrr::map(taxlist, 
            ~ mpa_string(unlist(.x), mpa_tib))) %>%
    unnest(tax_string) %>%
    dplyr::select(-taxlist) %>%
    dplyr::select(taxid, tidyselect::any_of(all_levels), tax_string)

# Check
stopifnot(nrow(taxtib) == nrow(summ_tb))
stopifnot(sum(summ_tb$taxid != taxtib$taxid)==0)
#left_join(summ_tb, taxtib)


# Select only subset of levels for the standardized
# output table
mpa_tib %<>%
    dplyr::filter(colname %in% colnames(taxtib))
nmap <- setNames(mpa_tib$new_colname, mpa_tib$colname)
sub_tax <- taxtib %>%
    dplyr::select(taxid, any_of(mpa_tib$colname), tax_string) 
sub_tax %<>%
    dplyr::rename_at(names(nmap), ~ nmap[.])

summ_tb %>%
    dplyr::select(taxname, taxid) %>%
    dplyr::left_join(sub_tax) %>%
    write_tsv(arguments$outfile_standard)
    
summ_tb %>%
    dplyr::select(taxname, taxid) %>%
    dplyr::left_join(taxtib) %>%
    readr::write_tsv(arguments$outfile_full)
