#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Define project
args = commandArgs(trailingOnly=TRUE)
proj = args[1]

# Define infiles
gene_file = "data/Ralstonia_genes.tab"
barc_file = paste("results/projects/", proj, "/", proj, ".poolcount", sep="")
meta_file = paste("data/projects/", proj, ".metadata.tab", sep="")

# Load data
gene = read_tsv(gene_file)
barc = read_tsv(barc_file)
meta = read_tsv(meta_file)

# Consider only barcodes inserted into the central 10% to 90% of the genes
gpos = gene %>%
  # Select only the necessary information from the genes file
  select(locusId, scaffoldId, begin, end) %>%
  apply(1, function(x){
    gene_length = as.numeric(x[4]) - as.numeric(x[3]) + 1
    # Consider only the central 10% to 90% of the genes
    gene_begin = round(as.numeric(x[3]) + 0.1*gene_length)
    gene_end = round(as.numeric(x[4]) - 0.1*gene_length)
    # Display all the positions inside the genes
    tibble(locusId = x[1], scaffold = x[2], pos = gene_begin:gene_end)
  }) %>%
  bind_rows()

# Associate barcodes with the genes according to their scaffolds and positions
barc = inner_join(barc, gpos) %>%
  # Gather the counts of each sample from a column to a row
  gather(
    ID, Counts, -locusId, -scaffold, -pos, -barcode, -rcbarcode, -strand
  ) %>%
  # Modify ID to match metadata file
  mutate(ID = unlist(lapply(str_split(ID, "_"), "[[", 2))) %>%
  # Add metadata for samples
  inner_join(meta)

# Select barcodes and genes with adequate coverage in time-zero samples
# 3 different dates, 2 replicates for each date, 6 time-zero samples in total
# Require:
# - at least 3 reads per strain
# - 30 reads per gene (consider only the adequate strains)
bar0 = barc %>%
  # Filter time-zero samples
  filter(startsWith(Condition, "Time0")) %>%
  # Sum per-strain counts across two replicate time-zero samples for each date
  group_by(barcode, Date) %>%
  mutate(n0 = sum(Counts)) %>%
  # At least 3 reads for each strain
  filter(n0 >= 3) %>%
  # At least 30 reads for each gene
  group_by(locusId, Date) %>%
  filter(sum(n0) >= 30)

# Filter to only the adequate barcodes in the barcode pool
barc = inner_join(barc, select(bar0, barcode, locusId, Date, n0) %>% distinct())

# Remove time-zero since it is now included as the variable n0
barc = filter(barc, !startsWith(Condition, "Time0"))

# Function to calculate strain variance
strain_variance = function(n, n0){(1/(1+n)+1/(1+n0))/(log(2)^2)}

# Calculate maximum weight (see its use below)
Weight_max = 1/strain_variance(20, 20)

# Strain fitness is the normalized log2 ratio of counts between the treatment
# sample and the reference time-zero sample.
# Gene fitness is the weighted average of the strain fitness.
gfit = barc %>%
  select(barcode, locusId, Date, Condition, Counts, Sample, n0) %>%
  # Calculate strain variance
  mutate(Strain_variance = strain_variance(Counts, n0)) %>%
  # The gene fitness is the weighted average of the strain fitness
  # Calculate the weight of each barcode
  # Set a ceiling on the weight, being that of a strain with 20 reads in each sample
  # The Weight = min(1/Strain_variance, Weight_max)
  mutate(
    Weight = ifelse(
      1/Strain_variance < Weight_max,
      1/Strain_variance,
      Weight_max
    )
  ) %>%
  group_by(locusId) %>%
  # Count the number of strains in each gene
  # (the number of barcodes inserted in each gene)
  mutate(Strains_per_gene = n_distinct(barcode)) %>%
  # Sum counts and n0 over all strains in all genes
  ungroup() %>%
  mutate(Sum_nafter = sum(Counts)) %>%
  mutate(Sum_n0 = sum(n0))

# Function to calculate strain fitness
strain_fitness = function(n, p, n0){log2(n + sqrt(p)) - log2(n0 + sqrt(1/p))}

# Calculate gene fitness for genes with different strain counts separately
# For genes with only 1 or 2 strains:
gLT3 = filter(gfit, Strains_per_gene < 3) %>%
  # The gene factor is one, and pseudocount will be Sum_nafter/Sum_n0
  mutate(Pseudocounts = 1 * Sum_nafter/Sum_n0) %>%
  # Calculate strain fitness
  mutate(Strain_fitness = strain_fitness(Counts, Pseudocounts, n0)) %>%
  # Calculate weighted gene fitness
  group_by(Sample, locusId) %>%
  mutate(Gene_fitness = sum(Weight * Strain_fitness)/sum(Weight))

# For genes with more than 2 strains:
gGE3 = filter(gfit, Strains_per_gene >= 3) %>%
  # Calculate a preliminary fitness per strain
  group_by(Sample, locusId) %>%
  mutate(Pre_strain_fitness = strain_fitness(Counts, 1, n0)) %>%
  # Calcalate preliminary fitness per gene as median of strains
  mutate(Pre_gene_fitness = median(Pre_strain_fitness))

# Calculate the median preliminary fitness across all genes and samples
pre_med = gGE3 %>%
  select(locusId, Sample, Pre_gene_fitness) %>%
  distinct() %>%
  pull(Pre_gene_fitness) %>%
  median()

gGE3 = gGE3  %>%
  # Normalize preliminary fitness per gene to a median of zero
  mutate(Pre_gene_fitness = Pre_gene_fitness - pre_med) %>%
  # Calculate pseudocounts using the gene factor 2 ^ Pre_gene_fitness
  mutate(Pseudocounts = (2^Pre_gene_fitness) * (Sum_nafter/Sum_n0)) %>%
  # Calculate strain fitness
  mutate(Strain_fitness = strain_fitness(Counts, Pseudocounts, n0)) %>%
  # Calculate weighted gene fitness
  mutate(Gene_fitness = sum(Weight * Strain_fitness)/sum(Weight))

# Concatenate the two tables of genes with different number of strains
gfit = bind_rows(
  select(gLT3, -Sum_nafter, -Sum_n0, -Pseudocounts),
  select(
    gGE3,
    -Sum_nafter, -Sum_n0, -Pseudocounts,
    -Pre_strain_fitness, -Pre_gene_fitness
  )
)

# Normalization
# Determine position in scaffold for each gene and save in Index column
gnrm = gene %>%
  # Calculate the Middle of the gene
  mutate(Middle = (begin + end)/2) %>%
  # Order genes by the Middle position
  arrange(Middle) %>%
  # For each scaffold...
  group_by(scaffoldId) %>%
  # ...add the Index (position) in that scaffold for every gene
  mutate(Index = 1:length(locusId)) %>%
  # Select only the necessary data for genes
  select(locusId, scaffoldId, Index)

# Join genes and fitness into one table
gnmf = inner_join(
  gnrm,
  # Create a table with the distinct gene fitness value for each sample and gene
  select(gfit, locusId, Sample, Gene_fitness) %>% distinct()
)

# Create a table with windows on which to calculate local medians
wind = gnrm %>%
  group_by(locusId, scaffoldId) %>%
  mutate(Window = list((Index-125):(Index+125))) %>%
  unnest(Window) %>%
  group_by(scaffoldId) %>%
  mutate(
    Window = case_when(
      Window < 1 ~ Window + max(Index),
      Window > max(Index) ~ Window - max(Index),
      Window >= 1 & Window <= max(Index) ~ Window
    )
  )

# Calculate the local median for each gene
locm = gnmf %>%
  select(-locusId) %>%
  rename(Window = Index) %>%
  inner_join(wind) %>%
  group_by(Sample, locusId) %>%
  summarise(LocalMedian = median(Gene_fitness))

# Calculate normalized gene fitness values
gnmf = gnmf %>%
  inner_join(locm) %>%
  mutate(Norm_fg_0 = Gene_fitness - LocalMedian)

# Define a function to find the mode in a vector using kernel density
get_mode = function(x){
  xd = density(x)
  return(xd$x[which.max(xd$y)])
}

# Subtract the mode for each scaffold
gnmf = gnmf %>%
  group_by(scaffoldId, Sample) %>%
  arrange(Index) %>%
  mutate(
    Mode = get_mode(Norm_fg_0),
    Norm_fg = Norm_fg_0 - Mode
  )

# Add the normalized gene fitness score to the original table
gfit = inner_join(gfit, select(gnmf, scaffoldId, locusId, Sample, Norm_fg))

# t-like test statistic
# Genes without at least 15 time-zero reads on each side are excluded
gmid = gene %>%
  select(locusId, scaffoldId, begin, end) %>%
  mutate(Middle = (end + begin)/2) %>%
  select(-begin, -end)

# Estimate a priori noise in gene fitness
side = inner_join(bar0, gmid) %>%
  # Determine which side of the middle each insertion site is located
  mutate(Side = ifelse(pos < Middle, "Left", "Right")) %>%
  # Ensure that at least 15 reads are found on the each side
  select(barcode, locusId, Date, Side, n0) %>%
  distinct() %>%
  group_by(locusId, Date, Side) %>%
  filter(sum(n0) >= 15) %>%
  ungroup() %>%
  select(-n0) %>%
  group_by(locusId, Date) %>%
  # Make sure that both sides are present (both side have at least 15 reads)
  mutate(Both = n_distinct(Side)) %>%
  filter(Both == 2) %>%
  select(-Both) %>%
  distinct() %>%
  # Add fitness data
  inner_join(gfit) %>%
  # Calculate the median fitness on each side
  group_by(locusId, Side) %>%
  summarise(Median_side = median(Strain_fitness)) %>%
  # Spread sides into two columns (Left and Right)
  spread(Side, Median_side) %>%
  # Calculate the absolute difference between the sides
  mutate(Abs_diff = abs(Right - Left))

# Calculate mad12, i.e. median absolute difference between the two halves
mad12 = median(side$Abs_diff)

# Calculate Vt (typical gene variance) using mad12
Vt = (mad12^2)/((2*0.674)^2)

# Calculate the naive strain variance Vn for all strains
gfit = gfit %>%
  mutate(Vn = strain_variance(Counts, n0))

# Calculate the median Vn across genes use to calculate Vt, for each sample
medv = gfit %>%
  # Use only genes used to calculate mad12 and Vt
  filter(locusId %in% side$locusId) %>%
  # Calculate median Vn for each sample
  group_by(Sample) %>%
  summarise(Median_Vn = median(Vn))

tfit = gfit %>%
  # Add median
  inner_join(medv) %>%
  # Calculate the prior estimate of gene variance Vg using median of Vn over
  # genes used to calculate Vt
  group_by(Sample, locusId) %>%
  mutate(Vg = Vt * Vn/Median_Vn) %>%
  # Calculate the weighted sum of squared differences between fi and fg
  mutate(Vi = Weight * (Strain_fitness - Norm_fg)^2) %>%
  # Calculate the estimated variance Ve
  mutate(Ve = (sum(Vi)/sum(Weight) + Vg) / Strains_per_gene) %>%
  # Calculate the t statistic
  mutate(t = Norm_fg/sqrt(0.1^2 + max(Ve, Vn))) %>%
  # If the absolute value of t is greater than 4 the fitness is significant
  mutate(Significant = as.numeric(abs(t) > 4))

# Save fitness table
write_tsv(
  tfit,
  gzfile(paste("results/projects/", proj, "/", proj, ".fitness.tab.gz", sep=""))
)
