#!/usr/bin/env Rscript

################################################################################
# post-process MAGeCK test output files
# part of zuberlab/mageck-nf pipeline at https://github.com/zuberlab/mageck-nf
#
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# started 2017/09/27
################################################################################

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# packages
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)

# command line arguments
args       <- commandArgs(trailingOnly = TRUE)
sgrna_file <- args[1]
gene_file  <- args[2]
ctrl_sgRNAs_file <- args[3]


# ------------------------------------------------------------------------------
# import
# ------------------------------------------------------------------------------
sgrna_raw <- read_tsv(sgrna_file, col_types = "ccccddddddddddc")
genes_raw <- read_tsv(gene_file, col_types = "cidddiiddddiid")
ctrl_sgRNAs_raw <- data.frame() %>%
  mutate(id = NA)

if(!is.na(ctrl_sgRNAs_file)){
  ctrl_sgRNAs_raw <- read_tsv(ctrl_sgRNAs_file, col_types = "c") %>%
    rename(id=1)
}


# ------------------------------------------------------------------------------
# format sgRNA results
# ------------------------------------------------------------------------------
sgrna <- sgrna_raw %>%
  rename_all(str_replace, "[.]", "_") %>%
  mutate_at(c("control_count", "treatment_count"), str_replace_all, "[/]", ";") %>%
  rename(id = sgrna,
         group = Gene,
         treatment_mean = treat_mean,
         adjusted_var = adj_var,
         lfc = LFC,
         fdr = FDR) %>%
  select(id,
         group,
         control_count,
         treatment_count,
         control_mean,
         treatment_mean,
         control_var,
         adjusted_var,
         lfc,
         score,
         p_low,
         p_high,
         p_twosided,
         fdr)

sgrna %>%
  write_tsv("guides_stats.txt")

# ------------------------------------------------------------------------------
# format gene results
# ------------------------------------------------------------------------------
genes <- genes_raw %>%
  rename_all(str_replace, "[|]", "_") %>%
  rename(group = id, guides = num) %>%
  gather(metric, value, neg_score:pos_lfc) %>%
  separate(metric, into = c("direction", "metric"), sep = "_") %>%
  spread(metric, value) %>%
  rename(p = `p-value`, guides_good = goodsgrna) %>%
  select(direction, group, guides, guides_good, lfc, score, p, fdr) %>%
  arrange(fdr)

genes %>%
  nest(-direction) %>%
  walk2(.x = .$data, .y = .$direction, .f = ~ write_tsv(.x, paste0("genes_", .y, "_stats.txt")))

# ------------------------------------------------------------------------------
# qc plots
# ------------------------------------------------------------------------------
### sgrna level power plot
qc1 <- sgrna %>%
  select(p_twosided, fdr) %>%
  arrange(p_twosided) %>%
  mutate(rank = seq(n())) %>%
  gather(metric, value, -rank) %>%
  ggplot(aes(x = rank, y = value, color = metric)) +
  geom_abline(aes(slope = 1 / max(rank), intercept = 0), linetype = 3) +
  geom_line() +
  theme_light() +
  ggtitle("Guide-level power curve")
ggsave(filename = "sgrna_power_curve.pdf", plot = qc1, device = "pdf", width = 6, height = 5)

### sgrna level MA plot
qc2 <- sgrna %>%
  select(id, control_mean, lfc, fdr) %>%
  mutate(fdr = if_else(fdr < 1e-6, 1e-6, fdr)) %>%
  ggplot(aes(x = control_mean, y = lfc)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(color = -log10(fdr)), size = 0.75) +
  scale_color_gradientn(colors = rainbow(7)) +
  geom_point(data = . %>% subset(id %in% ctrl_sgRNAs_raw$id),color = "black", size = 0.75) +
  scale_x_log10() +
  theme_light() +
  ggtitle("Guide-level MA plot")
ggsave(filename = "sgrna_ma_plot.pdf", plot = qc2, device = "pdf", width = 6, height = 5)

### gene level power curve
qc3 <- genes %>%
  select(direction, p, fdr) %>%
  mutate(direction = if_else(direction == "neg", "negative selection", "positive selection")) %>%
  arrange(p) %>%
  group_by(direction) %>%
  mutate(rank = seq(n())) %>%
  ungroup %>%
  gather(metric, value, p:fdr) %>%
  ggplot(aes(x = rank, y = value, color = metric)) +
  geom_abline(aes(slope = 1 / max(rank), intercept = 0), linetype = 3) +
  geom_line() +
  facet_wrap(~ direction) +
  theme_light() +
  ggtitle("Gene-level power curve")
ggsave(filename = "gene_power_curve.pdf", plot = qc3, device = "pdf", width = 9, height = 5)
