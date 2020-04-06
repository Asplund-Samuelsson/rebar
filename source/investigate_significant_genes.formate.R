#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Define project
args = commandArgs(trailingOnly=TRUE)
proj = args[1]

# Define infiles
tfit_file = paste(
  "results/projects/", proj, "/", proj, ".fitness.tab.gz", sep=""
)
anno_file = "data/Ralstonia_annotations.tab"
feat_file = "data/Ralstonia_features.tab"
rnks_file = "data/redMAGPIE_CBB_feature_ranks.tab"
subt_file = "data/redMAGPIE_Ralstonia_subtree_sign_feature_correlation.tab"
enri_file = "data/redMAGPIE_significantly_enriched_features.tab"

# Load data
tfit = read_tsv(tfit_file)
anno = read_tsv(anno_file)
feat = read_tsv(feat_file) %>% mutate(Feature = str_replace(Feature, ":", ""))
rnks = read_tsv(rnks_file, col_types=cols("c", "c", "i", "c", "d"))
subt = read_tsv(subt_file)
enri = read_tsv(enri_file)

# Reduce information
tfit = tfit %>%
  select(locusId, Date, Condition, Sample, Gene_fitness, t, Significant) %>%
  distinct()

# Select only significant formate data
tfit = filter(tfit, Condition == "Formate", Significant == 1)

# Add annotations
tfit = left_join(tfit, anno)

# Confirm that replicates agree
tfit %>%
  group_by(locusId) %>%
  summarise(
    Negative = sum(Gene_fitness > 0),
    Positive = sum(Gene_fitness < 0)
  ) %>%
  mutate(Disagreement = as.numeric(Negative > 0 & Positive > 0)) %>%
  filter(Disagreement != 0)

# Calculate average fitness across replicates
afit = tfit %>%
  group_by(locusId, protein_name, COG, COG_ID, COG_System, COG_Process) %>%
  summarise(Gene_fitness = mean(Gene_fitness)) %>%
  mutate(
    Label = paste(
      str_trunc(protein_name, 40, "right"), " [", locusId, "]",
      sep=""
    )
  )

# Check ranks
rfit = afit %>%
  ungroup() %>%
  select(locusId) %>%
  inner_join(feat) %>%
  inner_join(rnks) %>%
  arrange(Rank) %>%
  mutate(
    Method = recode(
      Ranking, 'ACE' = 'A', 'Enrichment' = 'E', 'Random forest' = 'R'
    )
  ) %>%
  select(-Ranking, -Diff) %>%
  spread(Method, Rank) %>%
  group_by(locusId) %>%
  summarise(
    A = mean(A, na.rm=T),
    E = mean(E, na.rm=T),
    R = mean(R, na.rm=T)
  ) %>%
  ungroup() %>%
  mutate(
    A = ifelse(is.na(A), max(A, na.rm=T)+1, A),
    E = ifelse(is.na(E), max(E, na.rm=T)+1, E),
    R = ifelse(is.na(R), max(R, na.rm=T)+1, R)
  ) %>%
  rowwise() %>%
  mutate(RankSum = sum(A, E, R)) %>%
  ungroup() %>%
  arrange(RankSum) %>%
  right_join(afit)

# Check correlations
cor(rfit$Gene_fitness, rfit$A, use="complete.obs", method = "spearman")
cor(rfit$Gene_fitness, rfit$E, use="complete.obs", method = "spearman")
cor(rfit$Gene_fitness, rfit$R, use="complete.obs", method = "spearman")
cor(rfit$Gene_fitness, rfit$RankSum, use="complete.obs", method = "spearman")

# Check correlation in subtree
sfit = afit %>% inner_join(feat) %>% inner_join(subt)
cor(sfit$Gene_fitness, sfit$r, use="complete.obs", method = "spearman")

# Check correlation versus significantly enriched features
efit = afit %>% inner_join(feat) %>% inner_join(enri)
cor(efit$Gene_fitness, efit$Diff, use="complete.obs", method = "spearman")
cor.test(efit$Gene_fitness, efit$Rank, use="complete.obs", method = "spearman")

# Perform tests between groups defined by formate fitness

# Load fitness again
dfit = read_tsv(tfit_file)
dfit = dfit %>%
  # Reduce information
  select(locusId, Condition, Sample, Gene_fitness, Significant) %>%
  distinct() %>%
  # Select formate data only
  filter(Condition == "Formate") %>%
  # Calculate average fitness across replicates
  group_by(locusId, Significant) %>%
  summarise(Gene_fitness = mean(Gene_fitness))

# Test differences between significant and non-significant depletion
wilcox.test(Diff ~ Significant, enri %>% inner_join(feat) %>% inner_join(dfit))
wilcox.test(r ~ Significant, subt %>% inner_join(feat) %>% inner_join(dfit))

# Summarise positive connection
chit = bind_rows(
    enri %>%
      inner_join(feat) %>%
      inner_join(dfit) %>%
      mutate(Connection = ifelse(Diff > 0, "+", "-")) %>%
      select(Significant, Connection) %>%
      mutate(Method = "Enrichment"),
    subt %>%
      inner_join(feat) %>%
      inner_join(dfit) %>%
      mutate(Connection = ifelse(r > 0, "+", "-")) %>%
      select(Significant, Connection) %>%
      mutate(Method = "ACE")
  ) %>%
  mutate(Significant = ifelse(Significant == 1, "Essential", "Irrelevant")) %>%
  group_by(Significant, Method, Connection) %>%
  summarise(Count = length(Connection)) %>%
  spread(Significant, Count)

# This is significant
chisq.test(select(ungroup(chit), Essential, Irrelevant))

# What is significant then?
chit = mutate(chit, Ratio = Essential / Irrelevant)

# Determine what genes are significant under non-formate conditions
nonf = read_tsv(tfit_file) %>%
  filter(Condition != "Formate", Significant == 1) %>%
  pull(locusId) %>%
  unique()

# Focus on formate unfit genes
ffit = filter(afit, !(locusId %in% nonf))
# 98.4% are specific to the formate condition!

## GO ENRICHMENT ###############################################################

library(topGO)

# Load GO terms
rean_file = "data/Ralstonia_H16_genome_re_annotation.csv"
rean = read_csv(rean_file)

# Create input table for topGO
godf = rean %>%
  dplyr::select(locus_tag, gene_ontology_IDs) %>%
  dplyr::rename(
    ClassicID = locus_tag, Gene.ontology.IDs = gene_ontology_IDs
  ) %>%
  inner_join(anno) %>%
  dplyr::select(locusId, Gene.ontology.IDs) %>%
  mutate(Formate = ifelse(locusId %in% afit$locusId, 1, 0)) %>%
  distinct() %>%
  as.data.frame()

#' Convenience wrapper to TopGO package (Rahnenfueher et al.)
#'
#' This function carries out a TopGO gene ontology enrichment on a data set
#' with custom protein/gene IDs and GO terms. The function takes as main
#' input a data frame with three specific columns, cluster numbers, Gene IDs,
#' and GO terms. Alternatively, these can also be supplied as three individual
#' lists.
#'
#' @param df an (optional) data.frame with the three columns specified below
#' @param GeneID (character) The column containing gene IDs, alternatively a vector
#' @param Gene.ontology.IDs (character) The column containing a list of GO terms for each gene,
#'   alternatively a vector with same order and length as 'GeneID'
#' @param cluster (numeric, factor, character) the column containing a grouping variable,
#'   alternatively a vector with same order and length as 'GeneID'
#' @param selected.cluster (character) the name of the group that is to be comapred to the background.
#'   Must be one of 'cluster'
#' @param topNodes (numeric) the max number of GO terms (nodes) to be returned by the function.
#'
#' @return a data.frame with TOpGO gene enrichment results
#'
#' @export
# ------------------------------------------------------------------------------

GetTopGO <- function(df = NULL, GeneID = NULL, Gene.ontology.IDs = NULL,
  cluster = NULL, selected.cluster, topNodes = 50
) {

  # prepare data structures
  # if a data.frame is passed as main data structure it must contain
  # three specific columns with cluster numbers, Gene IDs and GO terms.
  # GO terms is a character vector with go IDs separated by '; '
  #
  # if three separate lists are provided, they must be of same length
  # and order of geneIDs, cluster numbers and GO IDs must correspond
  # to each other

  # function to collect and prepare input data
  generate.input <- function(cluster, GeneID, Gene.ontology.IDs) {

    # test for duplicate IDs
    if (any(duplicated(GeneID)))
      stop("no duplicated GeneIDs allowed")

    # collect input
    genelist <- cluster; names(genelist) <- GeneID
    geneID2GO <- strsplit(Gene.ontology.IDs, "; ?")
    names(geneID2GO) <- GeneID

    # return two lists of genes and GO terms
    list(genelist = genelist, geneID2GO = geneID2GO)

  }

  if (class(df) == "data.frame" &
    all(c("cluster", "GeneID", "Gene.ontology.IDs") %in% colnames(df))) {

    input <- with(df,
      generate.input(cluster, GeneID, Gene.ontology.IDs)
    )

  } else if (!is.null(cluster) & !is.null(GeneID) & !is.null(Gene.ontology.IDs)) {
    input <- generate.input(cluster, GeneID, Gene.ontology.IDs)

  } else
    stop("no data provided or data not sufficiently formatted")


  # create topGO object
  topGOdata <- new("topGOdata",
    allGenes = input$genelist,
    annot = topGO::annFUN.gene2GO,
    gene2GO = input$geneID2GO,
    geneSel = function(x) x == selected.cluster,
    ontology = "BP")

  # We can use e.g. two types of test statistics: Fisher’s exact test which is based
  # on gene counts (e.g. genes from one cluster), and a Kolmogorov-Smirnov-
  # like test which computes enrichment based on gene scores (p-values).
  # Kolmogorov-Smirnov is only valid when p-values are provided, not a
  # set of interesting genes (e.g. a cluster)!
  resultFisherClassic <- topGO::runTest(topGOdata, algorithm = "classic", statistic = "fisher")
  resultFisherWeight <- topGO::runTest(topGOdata, algorithm = "weight", statistic = "fisher")
  resultFisherElim <- topGO::runTest(topGOdata, algorithm = "elim", statistic = "fisher")

  # collect all test results in one table
  GenTab <- topGO::GenTable(topGOdata,
    classicFisher = resultFisherClassic,
    weightedFisher = resultFisherWeight,
    elimFisher = resultFisherElim,
    orderBy = "elimFisher", ranksOf = "elimFisher",
    topNodes = topNodes)

  # add gene names that are contained in the respective cluster/GO term
  GenTab$SigGenes <- sapply(GenTab$GO.ID, function(term){
    paste(collapse = ",",
      topGO::genesInTerm(topGOdata, term)[[1]][topGO::scoresInTerm(topGOdata, term)[[1]] == selected.cluster]
    )
  })
  GenTab
}

# the GetTopGO function takes as input either
#   1. a data frame with one row per gene as only presence in cluster is
#      important here. Three columns are obligatory, "cluster", "GeneID", and
#      "Gene.ontology.IDs" with GO terms separated by '; '
#   2. alternatively three vectors or lists corresponding to
#      the three columns can be provided

TopGoResult = lapply(1, function(i) {
    with(godf,
      GetTopGO(df = NULL,
        cluster = as.numeric(Formate),
        GeneID = locusId,
        Gene.ontology.IDs,
        topNodes = 50,
        selected.cluster = i
      )
    )
  }) %>%

  # regarding p-value adjustment: the authors discourage from multiple
  # hypothesis testing and indeed, it turns most p-values insignificant
  setNames(1) %>%

  # turn into data frame with cluster as ID columns
  plyr::ldply(., .id = "Formate")

# save unfiltered TopGO result
write_csv(
  TopGoResult,
  paste(
    "intermediate/projects/", proj, "/", proj, ".unfiltered_TopGoResult.csv",
    sep=""
  )
)

# Get dispensability scores from REVIGO (http://revigo.irb.hr/)

# Load REVIGO results
rvgo_file = "intermediate/projects/rebar/rebar.REVIGO.csv"
rvgo = read_csv(rvgo_file)

# p-value has already been filtered by < 0.05 (manually)
gotb = rvgo %>%
  filter(dispensability <= 0.55) %>%
  arrange(`log10 p-value`) %>%
  dplyr::select(term_ID, description) %>%
  left_join(
    dplyr::select(
      TopGoResult, GO.ID, Annotated, Significant, Expected, elimFisher
    ) %>%
    dplyr::rename(term_ID = GO.ID)
  ) %>%
  filter(Annotated >= 5)

# Save results table
write_tsv(
  gotb,
  paste(
    "results/projects/", proj, "/", proj, ".formate_enriched_GOs.tab", sep=""
  )
)
