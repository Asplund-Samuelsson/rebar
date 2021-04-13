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
    strip.text = element_text(colour = "white", size = 10, margin = unit(rep(3,4), "points")),
    legend.text = element_text(colour = grey(0.4), size = 10),
    legend.title = element_text(colour = grey(0.4), size = 10),
  )
}

# QC PLOT 1: Number of reads per barcode/mutant, per sample
plot_read_count <- df_counts %>%
  ggplot(aes(x = log2(n_reads))) +
  geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  labs(x = expression("log"[2]*" reads per barcode")) +
  facet_wrap(~ sample) +
  custom_theme()

# QC PLOT 2: Number of barcodes/mutants per gene
plot_barcodes_gene <- fitness %>% ungroup %>%
  select(locusId, Strains_per_gene) %>%
  distinct %>% filter(Strains_per_gene < 40) %>%
  ggplot(aes(x = Strains_per_gene)) +
  geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  labs(x = "barcodes/mutants per gene") +
  custom_theme()

# QC PLOT 3: Average reads per gene (median of all samples)
plot_reads_gene <- fitness %>%
  group_by(locusId) %>%
  summarize(reads_per_gene_median = median(Counts)) %>%
  ggplot(aes(x = log2(reads_per_gene_median))) +
  geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  labs(x = expression("log"[2]*" reads per gene")) +
  custom_theme()


# PERFORM PCA
# ====================

# Reduce information
fitness = fitness %>%
  select(locusId, Date, Condition, ID, Norm_fg, t, Significant) %>%
  distinct() %>%
  mutate(ID = as.character(ID))

# Log-transform and center the data
wide = fitness %>%
  select(locusId, ID, Norm_fg) %>%
  spread(ID, Norm_fg) %>%
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
    select(fitness, locusId, Date, Condition, ID) %>% 
      distinct(), by = "ID")

# Calculate fraction of variance per PC
pcva = percent(fpca$sdev^2 / sum(fpca$sdev^2))[1:3]

plot_pca = ggplot(fplt, aes(x=PC1, y=PC2, label=ID, group=Condition, colour=Date)) +
  geom_line(colour="grey") +
  geom_point() +
  geom_text_repel(force=3, size=4) +
  labs(
    x=paste("PC1 (", pcva[1], ")", sep=""),
    y=paste("PC2 (", pcva[2], ")", sep="")
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