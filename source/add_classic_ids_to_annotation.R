#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Define infiles
anno_file = "data/Ralstonia_annotations.tab"
mich_file = "data/Ralstonia_H16_genome_re_annotation.csv"

# Load data
anno = read_tsv(anno_file)
mich = read_csv(mich_file)

mich = mich %>%
  select(locus_tag, new_locus_tag) %>%
  mutate(
    locusId = ifelse(
      new_locus_tag %in% anno$locusId,
      new_locus_tag, locus_tag
    )
  ) %>%
  rename(ClassicID = locus_tag) %>%
  select(-new_locus_tag) %>%
  distinct()

# Write updated annotation file
write_tsv(left_join(anno, mich), anno_file)
