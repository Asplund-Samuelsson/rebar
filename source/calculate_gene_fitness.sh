#!/usr/bin/env bash
#
# Script to merge sgRNA read counts from different samples into one table,
# perform DESeq2 analysis for differential abundance, and calculate fitness score

# Author: Michael Jahn
# Date: 2021-03-27

# optional input parameters
result=${result:-"./data/example/results/result.poolcount"}
poolfile=${poolfile:-"/ref/poolfile.tsv"}
gene_id=${gene_id:-"old_locus_tag"}
metadata=${metadata:-"./data/example/fastq/metadata.tsv"}
output_dir=${output_dir:-"./"}

# assign optional parameters that were passed with "--"
while [ $# -gt 0 ]; do
  if [[ $1 == *"--"* ]]; then
    param="${1/--/}"
    declare $param="$2"
  fi
  shift
done

# this bash script is just a wrapper around an R script that does all the work:
Rscript source/calculate_gene_fitness.R ${result} ${poolfile} \
  ${gene_id} ${metadata} ${output_dir}

# this R script performs PCA and generates summary plots
Rscript source/summary_plots.R ${result} ${output_dir}
