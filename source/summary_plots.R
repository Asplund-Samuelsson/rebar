#!/usr/bin/env Rscript
#
# Authors: 
# Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)
# Michael Jahn, KTH (michael.jahn@scilifelab.se)
# 
# Description: This script generates PCA and summary plots from BarSeq read counts

# LOAD PACKAGES
# ====================
#
cat("Loading required R packges: tidyverse, ggrepel, scales.\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(scales)
})

# load data from supplied arguments
args <- commandArgs(trailingOnly = TRUE)
df_counts <- read_tsv(args[1], col_types = cols())
load(paste0(args[2], "/fitness.Rdata"))

# reshape barcode counts to long format
df_counts <- df_counts %>% pivot_longer(
  cols = !all_of(c("barcode", "rcbarcode", "scaffold", "strand", "pos")), 
  names_to = "sample", values_to = "n_reads")


# SUMMARY PLOTS
# ====================

# define a custom ggplot2 theme (just for prettiness)
custom_colors = c("#E7298A", "#66A61E", "#E6AB02", "#7570B3", "#666666", "#1B9E77", "#D95F02", "#A6761D")
custom_theme <- function(base_size = 12, base_line_size = 1.0, base_rect_size = 1.0) {
  theme_light(base_size = base_size, base_line_size = base_line_size, base_rect_size = base_rect_size) + theme(
    plot.margin = unit(c(12,12,12,12), "points"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(colour = grey(0.4), linetype = "solid", lineend = "round"),
    axis.text.x = element_text(colour = grey(0.4), size = 10),
    axis.text.y = element_text(colour = grey(0.4), size = 10),
    panel.grid.major = element_line(size = 0.6, linetype = "solid", colour = grey(0.9)),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = grey(0.4), fill = NA, size = 1.0),
    panel.background = element_blank(),
    strip.background = element_rect(fill = grey(0.4), colour = grey(0.4)),
    strip.text = element_text(colour = "white", size = 10, margin = unit(rep(1,4), "points")),
    legend.text = element_text(colour = grey(0.4), size = 10),
    legend.title = element_text(colour = grey(0.4), size = 10),
  )
}

# QC PLOT: Total number of mapped reads per sample
plot_total_mapped_reads <- df_counts %>%
  group_by(sample) %>% summarize(n_reads = sum(n_reads)) %>%
  ggplot(aes(x = sample, y = n_reads)) +
  coord_flip() +
  geom_col(fill = custom_colors[1], alpha = 0.7) +
  labs(x = "", y = "total number of mapped reads") +
  custom_theme()

# QC PLOT: Number of reads per barcode/mutant, per sample
plot_read_count <- df_counts %>%
  ggplot(aes(x = log2(n_reads))) +
  geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  labs(x = expression("log"[2]*" reads per barcode")) +
  facet_wrap(~ sample) +
  custom_theme()

# QC PLOT: top 20 most abundant barcodes, per sample
plot_top_barcodes <- df_counts %>%
  group_by(sample) %>% arrange(sample, desc(n_reads)) %>% 
  mutate(rank = seq_along(barcode)) %>%
  filter(between(rank, 1, 10)) %>%
  ggplot(aes(x = factor(rank), y = n_reads)) +
  geom_col(fill = custom_colors[1], alpha = 0.7, width =1) +
  labs(y = "n reads", x = "barcode ranked by abundance") +
  facet_wrap(~ sample) +
  custom_theme()

# QC PLOT: Number of barcodes/mutants per gene
plot_barcodes_gene <- fitness %>% ungroup %>%
  select(locusId, strains_per_gene) %>%
  distinct %>% filter(strains_per_gene < 40) %>%
  ggplot(aes(x = strains_per_gene)) +
  geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  labs(x = "barcodes/mutants per gene") +
  custom_theme()

# QC PLOT: Average reads per gene (median of all samples)
plot_reads_gene <- fitness %>%
  group_by(locusId, barcode) %>%
  summarize(median_reads_per_bc = median(counts, na.rm = TRUE), .groups = "drop_last") %>%
  summarize(sum_reads_per_gene = sum(median_reads_per_bc, na.rm = TRUE)) %>%
  ggplot(aes(x = log2(sum_reads_per_gene))) +
  geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  labs(x = expression("log"[2]*" reads per gene (sum of all barcodes)")) +
  custom_theme()


# PERFORM PCA
# ====================

# Reduce information
fitness = fitness %>%
  select(locusId, date, condition, ID, norm_gene_fitness) %>%
  distinct() %>%
  mutate(ID = as.character(ID))

# Log-transform and center the data
wide = fitness %>% group_by(locusId) %>%
  select(locusId, ID, norm_gene_fitness) %>%
  spread(ID, norm_gene_fitness) %>%
  as.data.frame() %>%
  na.omit()

rownames(wide) = wide$locusId

fmat = select(wide, -locusId) %>%
  as.matrix() %>%
  t() %>%
  scale(center=T, scale=F) %>%
  t()

# Perform PCA
fpca = prcomp(fmat)

# Create plotting dataframes
fplt = as.data.frame(fpca$rotation)
fplt$ID = rownames(fplt)

# Add information about replicates, conditions, and dates
fplt = fplt %>%
  as_tibble() %>%
  select(ID, PC1, PC2) %>%
  inner_join(
    select(fitness, locusId, date, condition, ID) %>% 
      distinct(), by = "ID")

# Calculate fraction of variance per PC
pcva = percent(fpca$sdev^2 / sum(fpca$sdev^2))[1:3]


plot_pca = ggplot(select(fplt, -locusId) %>% distinct, 
  aes(x = PC1, y = PC2, label = ID, group = condition, colour = date)) +
  geom_line(colour = "grey") +
  geom_point() +
  geom_text_repel(force = 3, size = 4) +
  labs(
    x = paste("PC1 (", pcva[1], ")", sep = ""),
    y = paste("PC2 (", pcva[2], ")", sep = "")
  ) +
  custom_theme()


# EXPORT SUMMARY PLOTS
# ==============================

# export function
save_plots <- function(pl) {
  pdf(file = paste0(args[2], pl, ".pdf"), paper = "a4")
  print(get(pl))
  dev.off()
  png(filename = paste0(args[2], pl, ".png"), width = 1200, height = 1200, res = 120)
  print(get(pl))
  dev.off()
}

cat("Saving plots to", args[2], ".\n")
for (pl in grep(pattern = "^plot\\_", ls(), value = TRUE)) {
  save_plots(pl)
}
cat(" ---------------------------------\n", "Export of plots completed.\n")