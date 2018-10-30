#!/usr/bin/env Rscript

################################################################################
# Quantile normalize MAGeCK count files
# part of zuberlab/mageck-nf pipeline at https://github.com/zuberlab/mageck-nf
#
# Tobias Neumann
# Institute of Molecular Pathology (IMP), Vienna, Austria
# started 2018/10/30
################################################################################

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# packages
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(preprocessCore)

# command line arguments
args        <- commandArgs(trailingOnly = TRUE)
file_counts <- args[1]

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------
# prefilter_counts
prefilter_counts <- function(df, cols, min_count = 50) {
  
  if (!length(intersect(cols, names(df))) == length(cols)) {
    stop("`cols` not in column names of count matrix", call. = FALSE)
  }
  
  if (!all(sapply(df[, cols], is.numeric))) {
    stop("`cols` contains non-numeric columns", call. = FALSE)
  }
  
  df %>%
    dplyr::filter_at(vars(cols), all_vars(. >= min_count))
}

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------

counts = read_tsv(file_counts)
qn = normalize.quantiles(as.matrix(counts %>% select(-id, -group)),copy=TRUE)
out = data.frame(cbind(counts %>% select(id, group), qn))
names(out) = names(counts)

out %>%
	format_tsv %>%
	cat
  