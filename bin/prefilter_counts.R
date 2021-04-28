#!/usr/bin/env Rscript

################################################################################
# pre-process MAGeCK count files
# part of zuberlab/mageck-nf pipeline at https://github.com/zuberlab/mageck-nf
#
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# started 2018/01/05
################################################################################

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# packages
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(rlang)

# command line arguments
args        <- commandArgs(trailingOnly = TRUE)
file_counts <- args[1]
controls    <- args[2]
treatments  <- args[3]
min_count   <- args[4]

# format command line arguments
min_count <- as.integer(min_count)

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------
# prefilter_counts
prefilter_counts <- function(df, cols, min_count = 50) {
	
  if (length(grep(",",cols)) > 0 && length(grep("\\(",cols)) == 0) {
	  
	  cols  <- unlist(strsplit(cols, split = ","))
  
	  if (!length(intersect(cols, names(df))) == length(cols)) {
	    stop("`cols` not in column names of count matrix", call. = FALSE)
	  }
	  
	  if (!all(sapply(df[, cols], is.numeric))) {
	    stop("`cols` contains non-numeric columns", call. = FALSE)
	  }
	  
	  df %>%
	    dplyr::filter_at(vars(cols), all_vars(. >= min_count))

  } else {
	  
	  df %>%
	    dplyr::mutate(res := !!parse_quosure(cols)) %>%
        dplyr::filter(res >= min_count) %>% 
		dplyr::select(-res)
	  
  }
}

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
ctrl  <- unlist(strsplit(controls, split = ","))
trmnt  <- unlist(strsplit(treatments, split = ","))
readr::read_tsv(file_counts) %>%
  dplyr::select(c(1,2), ctrl, trmnt) %>%
  prefilter_counts(cols = controls, min_count = min_count) %>%
  tidyr::drop_na() %>%
  format_tsv %>%
  cat
  