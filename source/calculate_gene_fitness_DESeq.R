#!/usr/bin/env Rscript
#
# R script to summarize sequencing read counts obtained from a transposon library
#
# Author: Michael Jahn
# Description: The purpose of this script is to summarize count tables per sample
# in one main table, add statistical metrics for a pairwise sample comparison using DESeq2,
# and calculate fitness scores for each gene and condition.

# LOAD PACKAGES
# ====================
#
cat("Loading required R packges: DESeq2, DescTools, Hmisc, tidyverse.\n")
suppressPackageStartupMessages({
  library(DESeq2)
  library(DescTools)
  library(Hmisc)
  library(tidyverse)
})

# supplied input directories
args = commandArgs(trailingOnly = TRUE)
counts_dir <- args[1]
ref_dir <- args[2]
gene_col <- args[3]
metadata_dir <- args[4]
output_dir <- args[5]

# DATA PREPARATION
# ====================
#
# Step 1: Load sample layout sheet - file names must be row names
df_metadata <- read_tsv(metadata_dir, col_types = cols()) %>%
  mutate(file_name = gsub(".fastq.gz$", "", file_name)) %>%
  mutate(group = factor(`group`)) %>%
  column_to_rownames("file_name")
cat("Input:", nrow(df_metadata), "files listed in meta data table.\n")
stopifnot(is.numeric(df_metadata$time))

# Step 2: Load read counts
df_counts <- read_tsv(counts_dir, col_types = cols()) %>%
  pivot_longer(names_to = "file_name", values_to = "numreads", 
    cols = all_of(rownames(df_metadata))) %>%
  select(file_name, barcode, numreads)

# Step 3: Load gene reference
df_ref <- read_tsv(ref_dir, col_types = cols()) %>%
  rename(locusId = all_of(gene_col))

# print overview information to console
cat("Number of barcodes detected in n samples:\n")
df_counts %>% group_by(barcode) %>%
  summarize(barcodes_detected_in_samples = sum(numreads > 0)) %>%
  count(barcodes_detected_in_samples) %>%
  mutate(percent_total = n/sum(n)*100) %>%
  arrange(desc(barcodes_detected_in_samples)) %>%
  print

# DIFFERENTIAL ABUNDANCE
# ======================
#
# DESeq2 can be used to obtain fold changes and significance metrics
# for condition-wise comparisons, for details see publication:
# Love, M.I., Huber, W., Anders, S. Genome Biology, 15:550, 2014.
# (https://doi.org/10.1186/s13059-014-0550-8)
cat("Running DESeq2 for pairwise comparison.\nWarning: this step can be time and computation-intense.\n")

# 1. Read count matrix
# data frame must be reshaped to a 'counts matrix' with genes as rows
# and samples (conditions) as columns.
counts <- df_counts %>%
  # spread condition over columns and barcodes over rows
  pivot_wider(names_from = file_name, values_from = numreads) %>%
  # remove barcode column, replace NA with 0
  mutate_at(vars(-1), function(x) coalesce(x, 0)) %>%
  # add row_names from column 
  column_to_rownames("barcode")

# 2. Meta data
# Meta data is required to carry out the actual DESeq2 analysis
# by 'contrasting' (comparing) selected conditions to each other.
# We check that the order of file names corresponds to colnames of counts
stopifnot(colnames(counts) == row.names(df_metadata))

# 3. Perform DESeq2 analysis
# WARNING: This step is computation-intense and can take several hours
# for a large data set
DESeq_result <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = df_metadata,
  design = ~ group) %>%
  DESeq

# print overview about tested comparisons
# resultsNames(DESeq_result)

# The syntax to call DESeq2's `results(...)` function is to use one pair of 
# contrasts `contrast("variable", "level1", "level2")`. To automate this, 
# a list of condition and reference pairs is set up from meta data
combinations <- df_metadata %>%
  select(group, reference_group) %>%
  filter(group != reference_group) %>%
  mutate(across(.cols = everything(), .fns = as.character)) %>%
  distinct %>% as.list %>% transpose

# extract results for desired combinations
DESeq_result_table <- lapply(combinations, function(l) {
  results(DESeq_result, contrast = c("group", l$group, l$reference_group),
    parallel = TRUE, tidy = TRUE) %>% 
    as_tibble %>% mutate(group = l$group) %>% rename(barcode = row)
  }) %>% bind_rows

# MERGE DESEQ RESULTS
# ======================
#
# merge DESeq result table with meta data
DESeq_result_table <- select(df_metadata, -replicate, -ID) %>% 
  distinct %>% as_tibble %>%
  full_join(DESeq_result_table, by = "group") %>%
  
  # complete missing combinations of variables, here mostly all log2FC
  # values (0) for the reference conditions
  complete(barcode, nesting(condition, date, time, group, reference_group)) %>%
  mutate(
    log2FoldChange = replace_na(log2FoldChange, 0),
    lfcSE = replace_na(lfcSE, 0),
    pvalue = replace_na(pvalue, 1),
    padj = replace_na(padj, 1)
  )

# CALCULATE FITNESS SCORE
# =======================
#
# Here we define fitness score as the area under/over the curve for log2 fold change
# over time. Enrichment will result in a positive score, depletion 
# in a negative score. The fitness score is normalized to the maximum time 
# for a particular condition, and is therefore independent of the duration
# of the cultivations. Requires at least 2 time points
if (length(unique(DESeq_result_table$time)) > 1) {
  DESeq_result_table <- DESeq_result_table %>%
    arrange(barcode, condition, time) %>%
    group_by(barcode, condition) %>%
    mutate(fitness = DescTools::AUC(time, log2FoldChange)/max(time))
}

# SUMMARIZE TO GENE-WISE COMPARISON
# =================================
#
df_counts <- left_join(
    df_counts, 
    rownames_to_column(df_metadata, "file_name"),
    by = "file_name") %>% 
  group_by(barcode, condition, date, time, group) %>%
  summarize(counts = sum(numreads), .groups = "drop")
  
fitness <- DESeq_result_table %>%
  left_join(df_counts, by = c("barcode", "condition", "date", "time", "group")) %>%
  # merge with gene/barcode mapping reference
  left_join(select(df_ref, -rcbarcode, -n, -nTot), by = "barcode") %>%
  # remove non-central barcodes
  filter(central) %>%
  # add final number of quantified barcodes/strains per gene
  group_by(locusId) %>%
  mutate(strains_per_gene = length(unique(barcode))) %>%
  group_by(locusId, condition) %>%
  mutate(norm_gene_fitness = median(fitness, na.rm = TRUE)) %>%
  ungroup %>% mutate(significant = padj <= 0.01) %>%
  # rename some columns simply for compatibility with Wetmore et al. protocol
  # for example ID/group, and also fitness scores have slightly different meanings
  rename(ID = group, strain_fitness = fitness, t = stat, 
    log2FC = log2FoldChange) %>%
  select(
    barcode, locusId, scaffold, date, time, ID, condition, counts,
    strains_per_gene, strain_fitness, norm_gene_fitness, log2FC, lfcSE, t, padj, significant
  )

# Log2 FC, fitness score and p-values were calculated barcode-wise.
# Here we summarize some of these statistics as median over all barcodes/strains
# For the t and p-value we report min, max and median
fitness_gene <- fitness %>%
  group_by(locusId, scaffold, date, time, condition, strains_per_gene) %>%
  summarise(.groups = "drop",
    counts = sum(counts),
    log2FC = median(log2FC, na.rm = TRUE),
    norm_gene_fitness = unique(norm_gene_fitness),
    t_min = min(t),
    t_median = median(t),
    t_max = max(t),
    p_min = min(padj),
    p_median = median(padj),
    p_max = max(padj)
  )

# EXPORT PROCESSED DATA
# =====================
#
# Save result tables to output folder, in Rdata format
cat("Saving 'fitness.Rdata' and 'fitness_gene.Rdata' to", output_dir,
  ".\n ---------------------------------\n", "Fitness calculation completed.\n\n")
save(fitness, file = paste0(output_dir, "/fitness.Rdata"))
save(fitness_gene, file = paste0(output_dir, "/fitness_gene.Rdata"))