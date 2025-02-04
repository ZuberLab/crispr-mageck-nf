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
library(ggplot2)

# command line arguments
args        <- commandArgs(trailingOnly = TRUE)
file_counts <- args[1]
filter_columns    <- args[2]
min_count   <- args[3]
estimate_min_count_from_samples   <- args[4]
control_sample_id   <- args[5]
treat_sample_id   <- args[6]
variance_sample_id   <- args[7]
mode <- args[8]
norm_method <- args[9]


# format command line arguments
min_count <- as.integer(min_count)
estimate_min_count_from_samples <- as.integer(estimate_min_count_from_samples)
# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------
# prefilter_counts
prefilter_counts <- function(df, filter_cols, control_cols, treat_cols, min_count = 50, estimate_min_count_from_samples=1, mode="default") {
  
  if(estimate_min_count_from_samples==1){
    #calc min counts for filtering: 5% of median counts 
    min_count <- median(rowSums(df[,match(control_cols,names(df))], na.rm = TRUE)/length(control_cols)) * 0.05
    message(paste0("Filter for >= ", min_count, " counts in reference sample."))
    
  }
  
  if (length(grep(",",filter_cols)) > 0 && length(grep("\\(",filter_cols)) == 0) {
    filter_cols  <- unlist(strsplit(filter_cols, split = ","))
  
	  if (!length(intersect(filter_cols, names(df))) == length(filter_cols)) {
	    stop("`filter_cols` not in column names of count matrix", call. = FALSE)
	  }
	  
	  if (!all(sapply(df[, filter_cols], is.numeric))) {
	    stop("`filter_cols` contains non-numeric columns", call. = FALSE)
	  }
	  
    buff <- df %>%
	    dplyr::filter_at(vars(filter_cols), all_vars(. >= min_count))

  } else {

	  buff <- df %>%
	    dplyr::mutate(res := !!parse_quo(filter_cols, env = global_env())) %>%
      dplyr::filter(res >= min_count) %>%
	    dplyr::select(-res)
	  
  }
  
  if(mode=="enrichment"){
    buff <- buff %>%
      # Calculate average counts for treatment and control samples
      mutate(
        avg_treatment = rowMeans(select(., all_of(treat_cols)), na.rm = TRUE),
        avg_control = rowMeans(select(., all_of(control_cols)), na.rm = TRUE)
      ) %>%
      # Filter rows where the average treatment count >= average control count
      filter(avg_treatment >= avg_control) %>%
      # Remove the temporary columns
      select(-avg_treatment, -avg_control)
  }
  
  #remove rows where there is NA instead of count (avoid errors with mageck)
  buff[rowSums(is.na(buff))==0,]
  
}

# ------------------------------------------------------------------------------
# process
# ------------------------------------------------------------------------------
#get all needed columns
filter_cols_buff  <- stringr::str_replace(pattern = "pmax\\(|\\)", string =unlist(strsplit(filter_columns, split = ",")), replacement = "")
columns<-c("id", "group", unlist(filter_cols_buff),unlist(strsplit(control_sample_id, split = ",")),unlist(strsplit(treat_sample_id, split = ",")),unlist(strsplit(variance_sample_id, split = ",")))
#remove undefined columns
columns <- columns[!columns == "empty"] %>% unique
counts <- read_tsv(file_counts)
#select needed columns
counts<-counts[, columns]
#filter non-represented guides
data <- counts %>% dplyr::select(-id, -group)
index_represented <- rowMeans(data) > 0
counts <- counts[index_represented,]

cpm_table <- counts %>%
  # Select only sample columns for CPM calculation
  mutate(across(-c(id, group), 
                ~ (.x / sum(.x)) * 1e7, 
                .names = "{col}"))

#plot counts
long_data <- cpm_table %>%
  pivot_longer(cols = -c(id, group), 
               names_to = "Sample", 
               values_to = "CPM")

# Create the boxplot with ggplot2
boxplot <- ggplot(long_data, aes(x = Sample, y = CPM)) +
  geom_boxplot(outlier.color = "red", fill = "lightblue") + # Boxplot
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "darkgreen", 
               aes(group = Sample)) + # Add mean points
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "darkorange", 
               aes(group = Sample)) + # Add median points
  theme_minimal() + # Minimal theme
  labs(title = "Boxplot of sgRNA CPMs per Sample",
       x = "Samples", 
       y = "CPM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

# Print the plot
ggsave(filename = "boxplot_cpm.pdf", plot = boxplot, device = "pdf", width = 12, height = 10)

#filter
cpm_filtered <- cpm_table %>%
  prefilter_counts(filter_cols = filter_columns, 
                   control_cols=unlist(strsplit(control_sample_id, split = ",")), 
                   treat_cols = unlist(strsplit(treat_sample_id, split = ",")), 
                   min_count, 
                   estimate_min_count_from_samples,
                   mode)

if(norm_method == "none"){
  counts <- cpm_filtered

}else{
  counts <- counts %>%
    filter(id %in% cpm_filtered$id)
}

counts %>%
  format_tsv %>%
  cat
  