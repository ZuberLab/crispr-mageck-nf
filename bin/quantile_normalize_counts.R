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
# process
# ------------------------------------------------------------------------------

counts = read_tsv(file_counts)
qn = normalize.quantiles(as.matrix(counts %>% select(-id, -group)),copy=TRUE)
out = data.frame(cbind(counts %>% select(id, group), qn))
names(out) = names(counts)

out %>%
	format_tsv %>%
	cat
  