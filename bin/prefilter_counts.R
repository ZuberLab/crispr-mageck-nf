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
min_count   <- args[3]
estimate_min_count_from_samples   <- args[4]

# format command line arguments
min_count <- as.integer(min_count)
estimate_min_count_from_samples <- as.integer(estimate_min_count_from_samples)
# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------
# prefilter_counts
prefilter_counts <- function(df, cols, min_count = 50, estimate_min_count_from_samples=1) {
  if (length(grep(",",cols)) > 0 && length(grep("\\(",cols)) == 0) {
	  
	  cols  <- unlist(strsplit(cols, split = ","))
  
	  if (!length(intersect(cols, names(df))) == length(cols)) {
	    stop("`cols` not in column names of count matrix", call. = FALSE)
	  }
	  
	  if (!all(sapply(df[, cols], is.numeric))) {
	    stop("`cols` contains non-numeric columns", call. = FALSE)
	  }
	  
	  if(estimate_min_count_from_samples==1){
	    #calc min counts for filtering: 5% of median counts 
	    
	    #mean version: 
	    #min_count <- round((sum(df[,match(cols,names(df))], na.rm = TRUE)/length(cols) /nrow(df)) * 0.05)
	    
	    #median
	    min_count <- round((median(rowSums(df[,match(cols,names(df))], na.rm = TRUE)/length(cols)) * 0.05))
	    message(paste0("Filter for >= ", min_count, " counts in reference sample."))
	    
	  }
	  
	  df %>%
	    dplyr::filter_at(vars(cols), all_vars(. >= min_count))

  } else {
    
    if(estimate_min_count_from_samples==1){
      cols_buff  <- stringr::str_replace(pattern = "pmax\\(|\\)", string =unlist(strsplit(cols, split = ",")), replacement = "")
      #calc min counts for filtering: 5% of median counts 
      
      #mean version: 
      #min_count <- round((sum(df[,match(cols_buff,names(df))], na.rm = TRUE)/length(cols_buff) /nrow(df)) * 0.05)
      
      #median
      min_count <- round((median(rowSums(df[,match(cols_buff,names(df))], na.rm = TRUE)/length(cols_buff)) * 0.05))
      message(paste0("Filter for >= ", min_count, " counts in reference sample."))
      
    }
	  
	  buff <- df %>%
	    dplyr::filter_at(vars(cols), all_vars(. >= min_count))
# 	    dplyr::mutate(res := !!parse_quosure(cols)) %>%
#       dplyr::filter(res >= min_count) %>% 
# 	    dplyr::select(-res)
	  
	  buff[,colSums(is.na(buff))==0]
	  
  }
}

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
read_tsv(file_counts) %>%
  prefilter_counts(cols = controls, min_count = min_count, estimate_min_count_from_samples=estimate_min_count_from_samples) %>%
  format_tsv %>%
  cat
  