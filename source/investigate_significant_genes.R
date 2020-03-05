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

# Load data
tfit = read_tsv(tfit_file)
anno = read_tsv(anno_file)

# Reduce information
tfit = tfit %>%
  select(locusId, Date, Condition, Sample, Gene_fitness, t, Significant) %>%
  distinct()

# Remove outliers and select only significant genes, then modify Isobutanol
tfit = tfit %>%
  filter(!(Sample %in% c("3a","3b")), Significant == 1) %>%
  mutate(
    Condition = ifelse(
      startsWith(Condition, "Isobutanol"), "Isobutanol", Condition
    )
  )

# Add annotations
tfit = left_join(tfit, anno)

# Confirm that replicates agree
tfit %>%
  group_by(Condition, locusId) %>%
  summarise(
    Negative = sum(Gene_fitness > 0),
    Positive = sum(Gene_fitness < 0)
  ) %>%
  mutate(Disagreement = as.numeric(Negative > 0 & Positive > 0)) %>%
  filter(Disagreement != 0)

# Calculate average fitness across replicates
afit = tfit %>%
  group_by(
    locusId, Condition, protein_name,
    COG, COG_ID, COG_System, COG_Process
  ) %>%
  summarise(Gene_fitness = mean(Gene_fitness)) %>%
  mutate(
    Label = paste(
      str_trunc(protein_name, 40, "right"), " [", locusId, "]",
      sep=""
    )
  ) %>%
  mutate(
    Label = factor(
      Label,
      levels = group_by(., Label) %>%
        summarise(f = sum(Gene_fitness)) %>%
        arrange(f) %>%
        pull(Label)
    )
  )

# Cluster conditions to make a nice order
cfit = afit %>%
  ungroup() %>%
  select(locusId, Condition, Gene_fitness) %>%
  spread(Condition, Gene_fitness) %>%
  as.data.frame()

rownames(cfit) = cfit$locusId
cfit = select(cfit, -locusId)

cfit[is.na(cfit)] = 0

afit = afit %>%
  ungroup() %>%
  mutate(
    Condition = factor(
      Condition,
      levels = colnames(cfit)[hclust(dist(t(cfit)))$order]
    )
  )

# Make a plot
gp = ggplot(afit, aes(y = Label, x = Condition, fill = Gene_fitness))
gp = gp + geom_tile()
gp = gp + scale_fill_gradient2(
  low="#af8dc3", mid="#f7f7f7", high="#7fbf7b", midpoint=0
)
gp = gp + facet_grid(
  COG_System + COG_Process ~ .,
  scales="free_y", space="free"
)
gp = gp + theme_bw()
gp = gp + theme(
  strip.background = element_blank(),
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
  strip.text.y = element_text(angle=0, hjust=0),
  axis.title = element_blank(),
  legend.position = c(-0.5,-0.04)
)

outfile = paste(
  "results/projects/", proj, "/", proj, ".sign_heatmap.pdf", sep=""
)

ggsave(outfile, gp, height=40, width=40, unit="cm")
