#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)
library(ggrepel)
library(scales)

# Define project
args = commandArgs(trailingOnly=TRUE)
proj = args[1]

# Define infiles
tfit_file = paste(
  "results/projects/", proj, "/", proj, ".fitness.tab.gz", sep=""
)

# Load data
tfit = read_tsv(tfit_file)

# Reduce information
tfit = tfit %>%
  select(locusId, Date, Condition, Sample, Gene_fitness, t, Significant) %>%
  distinct()

# Log-transform and center the data
wide = tfit %>%
  select(locusId, Sample, Gene_fitness) %>%
  spread(Sample, Gene_fitness) %>%
  as.data.frame() %>%
  na.omit()

rownames(wide) = wide$locusId

# Omit outliers
wide = select(wide, -`3a`, -`3b`)

fmat = select(wide, -locusId) %>%
  as.matrix() %>%
  t() %>%
  scale(center=T, scale=F) %>%
  t()

# Perform PCA
fpca = prcomp(fmat)

# Create plotting dataframes
fplt = as.data.frame(fpca$rotation)
fplt$Sample = rownames(fplt)

# Add information about replicates, conditions, and dates
fplt = fplt %>%
  as_tibble() %>%
  select(Sample, PC1, PC2) %>%
  inner_join(select(tfit, Date, Condition, Sample) %>% distinct()) %>%
  mutate(
    Condition = ifelse(
      startsWith(Condition, "Isobutanol"), "Isobutanol", Condition
    ),
    Label = ifelse(
      grepl("a", Sample), "", Condition
    )
  )

# Calculate fraction of variance per PC
pcva = percent(fpca$sdev^2 / sum(fpca$sdev^2))[1:3]

gp = ggplot(fplt, aes(x=PC1, y=PC2, label=Label, group=Condition, colour=Date))
gp = gp + geom_line(colour="grey")
gp = gp + geom_point()
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + scale_colour_manual(values=c("#f1a340", "#998ec3"))
gp = gp + labs(
            x=paste("PC1 (", pcva[1], ")", sep=""),
            y=paste("PC2 (", pcva[2], ")", sep="")
          )
gp = gp + theme_bw()

outfile = paste("results/projects/", proj, "/", proj, ".PCA.pdf", sep="")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)
