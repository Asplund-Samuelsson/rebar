library(tidyverse)

#define infiles
genes = "Ralstonia_genes.tab"
barcode_qi = "Barseq3.poolcount"
barcode_kyle = "Barseq2.poolcount"

#read infiles
genes = read_tsv(genes)
barcode_qi = read_tsv(barcode_qi)
barcode_kyle = read_tsv(barcode_kyle)

#only barcodes inserted into the central 10% to 90% of the genes are considered
#display all the gene positions
genes_position = genes %>%
  #select only the necessary information from the genes file
  select(locusId, scaffoldId, begin, end) %>%
  apply(1, function(x){
    gene_length = as.numeric(x[4]) - as.numeric(x[3]) + 1
    #consider only the central 10% to 90% of the genes
    gene_begin = round(as.numeric(x[3]) + 0.1*gene_length)
    gene_end = round(as.numeric(x[4]) - 0.1*gene_length)
    #display all the positions inside the genes
    tibble(locusId = x[1], scaffold = x[2], pos = gene_begin:gene_end)
  }) %>%
  bind_rows()


# A. For qi's data:
# (a) Select a subset of strains and genes that have adequate coverage in the time-zero samples.
#associate barcodes with the genes according to their scaffolds and positions
barcode_qi = inner_join(barcode_qi, genes_position) %>%
  #gather the counts of each sample from a column to a row
  gather(id, counts, -locusId, -scaffold, -pos, -barcode, -rcbarcode, -strand) %>%
  #add date
  mutate(date = recode(id, "IT001"="8/14","IT002"="8/14","IT003"="8/28","IT004"="8/14",
                       "IT005"="8/14","IT006"="8/28","IT007"="8/14","IT008"="8/14",
                       "IT009"="8/16","IT010"="8/14","IT011"="8/16","IT012"="8/14",
                       "IT013"="8/14","IT014"="8/14","IT015"="8/16","IT016"="8/16",
                       "IT017"="8/16","IT018"="8/16","IT019"="8/14","IT020"="8/16",
                       "IT021"="8/14","IT022"="8/14","IT023"="8/14","IT024"="8/14",
                       "IT025"="8/16","IT026"="8/16","IT027"="8/16","IT028"="8/16",
                       "IT029"="8/16","IT030"="8/16","IT031"="8/16","IT032"="8/16",
                       "IT033"="8/16","IT034"="8/16","IT035"="8/16","IT036"="8/16",
                       "IT037"="8/14","IT038"="8/14","IT039"="8/16","IT040"="8/16",
                       "IT041"="8/28","IT042"="8/28")) %>%
  #add condition
  #IT033,IT034(39a,39b) are from the same tube; IT035,IT036(40a,40b) are from the same tube;
  mutate(condition = recode(id, "IT001"="Fructose","IT002"="Succinate",
                            "IT003"="Formate","IT004"="Fructose",
                            "IT005"="Succinate","IT006"="Formate",
                            "IT007"="4-hydroxybenzoate","IT008"="4-hydroxybenzoate",
                            "IT009"="Propionic Acid","IT010"="Octanoic Acid",
                            "IT011"="Propionic Acid","IT012"="Octanoic Acid",
                            "IT013"="Oleic Acid","IT014"="Oleic Acid",
                            "IT015"="FructoseNaCl","IT016"="FructoseNaCl",
                            "IT017"="Glycerol","IT018"="Glycerol",
                            "IT019"="Isobutanol1","IT020"="Isobutanol2",
                            "IT021"="n-butanol","IT022"="n-butanol",
                            "IT023"="Ethanol","IT024"="Ethanol",
                            "IT025"="FructoseNaNO3","IT026"="FructoseNaNO3",
                            "IT027"="SuccinateNaNO3","IT028"="SuccinateNaNO3",
                            "IT029"="Octanoic Acid NaNO3","IT030"="Octanoic Acid NaNO3",
                            "IT031"="FructoseDMSO","IT032"="FructoseDMSO",
                            "IT033"="Fructose34C","IT034"="Fructose34C",
                            "IT035"="Fructose40C","IT036"="Fructose40C",
                            "IT037"="Time0_8/14","IT038"="Time0_8/14",
                            "IT039"="Time0_8/16","IT040"="Time0_8/16",
                            "IT041"="Time0_8/28","IT042"="Time0_8/28")) %>%
  #add sample name
  mutate(sample = recode(id, "IT001"="1a","IT002"="2a","IT003"="3a","IT004"="1b",
                         "IT005"="2b","IT006"="3b","IT007"="7a","IT008"="7b",
                         "IT009"="11b","IT010"="13a","IT011"="11a","IT012"="13b",
                         "IT013"="14a","IT014"="14b","IT015"="15a","IT016"="15b",
                         "IT017"="18a","IT018"="18b","IT019"="25a","IT020"="25b",
                         "IT021"="26a","IT022"="26b","IT023"="29a","IT024"="29b",
                         "IT025"="30a","IT026"="30b","IT027"="32a","IT028"="32b",
                         "IT029"="35a","IT030"="35b","IT031"="36a","IT032"="36b",
                         "IT033"="39a","IT034"="39b","IT035"="40a","IT036"="40b",
                         "IT037"="t0a","IT038"="t0b","IT039"="t02a","IT040"="t02b",
                         "IT041"="t03a","IT042"="t03b"))

#select barcodes and genes with adequate coverage in time-zero samples
#3 different dates, 2 replicates for each date, 6 time-zero samples in total
#at least 3 reads per strain, 30 reads per gene (consider only the adequate strians)
barcode_t0_qi = barcode_qi %>%
  #filter time-zero samples
  filter(id %in% c("IT037", "IT038", "IT039", "IT040", "IT041", "IT042")) %>%
  #sum the per-strain counts across two replicate time-zero samples on each date
  group_by(barcode, date) %>%
  mutate(n0 = sum(counts)) %>%
  #at least 3 reads for each strian
  filter(n0 >= 3) %>%
  #at least 30 reads for each gene
  group_by(locusId, date) %>%
  filter(sum(n0) >= 30)

#filter only the adequate barcodes in the barcode pool
barcode_qi = inner_join(barcode_qi,
                        select(barcode_t0_qi, barcode, locusId, date, n0) %>% distinct())

# (b) Strain fitness & Gene fitness
# Strain fitness is the normalized log2 ratio of counts between the treatment sample and the reference time-zero sample.
# Gene fitness is the weighted average of the strain fitness.
gene_fitness_qi = barcode_qi %>%
  select(barcode, locusId, date, condition, counts, sample, n0) %>%
  group_by(sample) %>%
  #the fraction of each barcode in each sample
  mutate(fraction = counts/sum(counts)) %>%
  #log2 fold change of each barcode in each sample before and after treatment
  mutate(log2FC = log2(counts/n0)) %>%
  #strain varience
  mutate(strain_variance = (1/(1+counts)+1/(1+n0))/log(4, base = exp(1)) ) %>%
  #the gene fitness is the weighted average of the strain fitness
  #calculate the weight of each barcode
  mutate(weight_cal = 1/strain_variance) %>%
  #set a ceiling on the weight, being that of a strain with 20 reads in each sample
  mutate(weight_max = 1/(1/(1+20)+1/(1+n0))/log(4, base = exp(1))) %>%
  #the weight = min (weight_cal, weight_max)
  mutate(weight = ifelse(weight_cal < weight_max, weight_cal, weight_max)) %>%
  group_by(locusId) %>%
  #the number of strains in each gene (the number of barcodes inserted in each gene)
  mutate(strains_per_gene = n_distinct(barcode)) %>%
  ungroup() %>%
  #the sum of counts and n0 over all strains in all genes
  mutate(sum_nafter = sum(counts)) %>%
  mutate(sum_n0 = sum(n0))

#calculate gene fitness for genes with different strains separately
#for genes with only 1 or 2 strains:
genes_1or2_strains_qi = filter(gene_fitness_qi, strains_per_gene < 3) %>%
  mutate(pseudocounts = sum_nafter/sum_n0) %>%
  mutate(strain_fitness = log2(counts + sqrt(pseudocounts)) - log2(n0 + sqrt(1/pseudocounts))) %>%
  group_by(sample, locusId) %>%
  mutate(gene_fitness = sum(weight * strain_fitness)/sum(weight))

#for genes with more than 2 strains:
genes_3_strains_qi = filter(gene_fitness_qi, strains_per_gene >= 3) %>%
  group_by(sample, locusId) %>%
  mutate(pre_strain_fitness = log2(counts + 1) - log2(n0 + 1)) %>%
  mutate(pre_gene_fitness = median(pre_strain_fitness)) %>%
  mutate(pseudocounts = (2^pre_gene_fitness) * (sum_nafter/sum_n0)) %>%
  mutate(strain_fitness = log2(counts + sqrt(pseudocounts)) - log2(n0 + sqrt(1/pseudocounts))) %>%
  mutate(gene_fitness = sum(weight * strain_fitness)/sum(weight))

#concatenate the two tables of genes with different number of strains:
gene_fitness_qi = bind_rows(select(genes_1or2_strains_qi, -weight_cal, -weight_max, -sum_nafter, -sum_n0, -pseudocounts),
                            select(genes_3_strains_qi, -weight_cal, -weight_max, -sum_nafter, -sum_n0, -pseudocounts, -pre_strain_fitness, -pre_gene_fitness))%>%
  mutate(identifier = "qi")

# (c) Normalization
#determine position in scaffold for each gene and save in Index column
genes_normalization = genes %>%
  #calculate the Middle of the gene
  mutate(Middle = (begin + end)/2) %>%
  #order genes by the Middle position
  arrange(Middle) %>%
  #for each scaffold...
  group_by(scaffoldId) %>%
  #...add the Index (position) in that scaffold for every gene
  mutate(Index = 1:length(locusId)) %>%
  #select only the necessary data for genes
  select(locusId, scaffoldId, Index)

#join genes and fitness into one table
normalization_qi = inner_join(genes_normalization,
                              #create a table with the distinct gene fitness value for each sample and gene
                              select(gene_fitness_qi, locusId, sample, gene_fitness) %>% distinct())

#create a table with windows on which to calculate local medians
windows = genes_normalization %>%
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

#calculate the local median for each gene
local_medians = normalization_qi %>%
  select(-locusId) %>%
  rename(Window = Index) %>%
  inner_join(windows) %>%
  group_by(sample, locusId) %>%
  summarise(LocalMedian = median(gene_fitness))

#calculate normalized gene fitness values
normalization_qi = normalization_qi %>%
  inner_join(local_medians) %>%
  mutate(norm_fg_0 = gene_fitness - LocalMedian)

#define a function to find the mode in a vector using kernel density
get_mode = function(x){
  xd = density(x)
  return(xd$x[which.max(xd$y)])
}

#subtract the mode for each scaffold
normalization_qi = normalization_qi %>%
  group_by(scaffoldId, sample) %>%
  arrange(Index) %>%
  mutate(
    Mode = get_mode(norm_fg_0),
    norm_fg = norm_fg_0 - Mode
  )

#add the normalized gene fitness score to the original table
gene_fitness_qi = inner_join(
  gene_fitness_qi, select(normalization_qi, scaffoldId, locusId, sample, norm_fg)
)

# (d) T-test
#genes without at least 15 time-zero reads on each side are excluded
genes_middle = genes %>%
  select(locusId, scaffoldId, begin, end) %>%
  mutate(middle = (end + begin)/2) %>%
  select(-begin, -end)

barcode_t0_side_qi = inner_join(barcode_t0_qi, genes_middle) %>%
  mutate(side = ifelse(pos < middle, "left", "right")) %>%
  group_by(locusId, date, side) %>%
  filter(sum(n0) >= 15) %>%
  ungroup() %>%
  select(barcode, locusId, date, side) %>%
  distinct() %>%
  group_by(locusId, date) %>%
  mutate(both = n_distinct(side)) %>%
  filter(both == 2) %>%
  select(-both) %>%
  distinct()

t_test_qi = inner_join(gene_fitness_qi, barcode_t0_side_qi) %>%
  group_by(locusId, sample, side) %>%
  mutate(median_side = median(strain_fitness)) %>%
  group_by(sample, locusId) %>%
  mutate(mad12 = diff(range(median_side))) %>%
  mutate(Vt = (mad12^2)/((2*0.674)^2)) %>%
  mutate(Vn = (1/(1+counts) + 1/(1+n0)) / log(4, base = exp(1))) %>%
  group_by(sample, locusId) %>%
  mutate(Vg = Vt * (Vn/median(Vn))) %>%
  mutate(Vi = weight * (strain_fitness - norm_fg)^2) %>%
  mutate(Ve = ((sum(Vi)/sum(weight)) + Vg)/ strains_per_gene) %>%
  mutate(t = norm_fg/sqrt(0.1^2 + max(Ve, Vn))) #%>%
  filter(abs(t)>4)

# B. For kyle's data:
# (a) Select a subset of strains and genes that have adequate coverage in the time-zero samples.
#associate barcodes with the genes according to their scaffold and positions
barcode_kyle = inner_join(barcode_kyle, genes_position) %>%
  #gather the counts of each sample from a column to a row
  gather(id, counts, -locusId, -scaffold, -pos, -barcode, -rcbarcode, -strand) %>%
  #add date
  mutate(date = recode(id, "IT001"="8/14","IT002"="8/14","IT003"="8/14","IT004"="8/14",
                       "IT005"="8/14","IT006"="8/14","IT007"="8/14","IT008"="8/14",
                       "IT009"="8/14","IT010"="8/14","IT011"="8/14","IT012"="8/14",
                       "IT013"="8/16","IT014"="8/16","IT015"="8/16","IT016"="8/14",
                       "IT017"="8/16","IT018"="8/16","IT019"="8/16","IT020"="8/16",
                       "IT021"="8/16","IT022"="8/16","IT023"="8/16","IT024"="8/16",
                       "IT025"="8/16","IT026"="8/16","IT027"="8/16","IT028"="8/16",
                       "IT029"="8/16","IT030"="8/16","IT031"="8/14","IT032"="8/14",
                       "IT033"="8/16","IT034"="8/16","IT035"="8/16","IT036"="8/16",
                       "IT037"="8/28","IT038"="8/28","IT039"="8/28","IT040"="8/28",
                       "IT041"="8/28","IT042"="8/28","IT043"="8/28","IT044"="8/28",
                       "IT045"="8/28","IT046"="8/28","IT047"="8/28")) %>%
  #add conditon
  mutate(condition = recode(id, "IT001"="Time0_8/14","IT002"="Time0_8/14",
                            "IT003"="Fructose","IT004"="Fructose",
                            "IT005"="Succinate","IT006"="Succinate",
                            "IT007"="4-hydroxybenzoate","IT008"="4-hydroxybenzoate",
                            "IT009"="Octanoic Acid","IT010"="Octanoic Acid",
                            "IT011"="Oleic Acid","IT012"="Oleic Acid",
                            "IT013"="Glycerol","IT014"="Glycerol",
                            "IT015"="Isobutanol2","IT016"="Isobutanol1",
                            "IT017"="n-butanol","IT018"="Prenol",
                            "IT019"="Time0_8/16","IT020"="Time0_8/16",
                            "IT021"="FructoseNaNO3","IT022"="FructoseNaNO3",
                            "IT023"="SuccinateNaNO3","IT024"="SuccinateNaNO3",
                            "IT025"="Octanoic Acid NaNO3","IT026"="Octanoic Acid NaNO3",
                            "IT027"="FructoseDMSO","IT028"="FructoseDMSO",
                            "IT029"="Malate","IT030"="Malate",
                            "IT031"="Ethanol","IT032"="Ethanol",
                            "IT033"="Propionic Acid","IT034"="Propionic Acid",
                            "IT035"="Butyric Acid","IT036"="Butyric Acid",
                            "IT037"="Time0_8/28","IT038"="Time0_8/28",
                            "IT039"="Formate","IT040"="Formate",
                            "IT041"="Malonate","IT042"="Malonate",
                            "IT043"="Palmitic Acid 1g/L + 10g/L DMSO","IT044"="Palmitic Acid 1g/L + 10g/L DMSO",
                            "IT045"="FACS Nile Red PHB positive population","IT046"="FACS Nile Red PHB middle population",
                            "IT047"="FACS Nile Red PHB-")) %>%
  #add sample name
  mutate(sample = recode(id, "IT001"="t0a","IT002"="t0b","IT003"="1a","IT004"="1b",
                         "IT005"="2a","IT006"="2b","IT007"="7a","IT008"="7b",
                         "IT009"="13a","IT010"="13b","IT011"="14a","IT012"="14b",
                         "IT013"="18a","IT014"="18b","IT015"="25b","IT016"="25a",
                         "IT017"="17a","IT018"="28a","IT019"="t02a","IT020"="t02b",
                         "IT021"="30a","IT022"="30b","IT023"="32a","IT024"="32b",
                         "IT025"="35a","IT026"="35b","IT027"="36a","IT028"="36b",
                         "IT029"="5a","IT030"="5b","IT031"="29a","IT032"="29b",
                         "IT033"="11a","IT034"="11b","IT035"="12a","IT036"="12b",
                         "IT037"="t03a","IT038"="t03b","IT039"="3a","IT040"="3b",
                         "IT041"="6a","IT042"="6b","IT043"="21a","IT044"="21b",
                         "IT045"="44a","IT046"="43a","IT047"="42a"))

#select barcodes and genes with adequate coverage in time-zero samples
#3 different date, 2 replicates for each date, 6 time-zero samples in total
#at least 3 reads per strain, 30 reads per gene (consider only the adequate strians)
barcode_t0_k = barcode_kyle %>%
  #filter time-zero samples
  filter(id %in% c("IT001", "IT002", "IT019", "IT020", "IT037", "IT038")) %>%
  #sum the per-strain counts across two replicate time-zero samples on each date
  group_by(barcode, date) %>%
  mutate(n0 = sum(counts)) %>%
  #at least 3 reads for each strian
  filter(n0 >= 3) %>%
  #at least 30 reads for each gene
  group_by(locusId, date) %>%
  filter(sum(n0) >= 30)

#filter only the adequate barcodes in the barcode pool
barcode_kyle = inner_join(barcode_kyle,
                          select(barcode_t0_k, barcode, locusId, date, n0) %>% distinct())

# (b) Strain fitness & Gene fitness
# Strain fitness is the normalized log2 ratio of counts between the treatment sample and the reference time-zero sample.
# Gene fitness is the weighted average of the strain fitness.
gene_fitness_k = barcode_kyle %>%
  select(barcode, locusId, date, condition, counts, sample, n0) %>%
  group_by(sample) %>%
  #the fraction of each barcodes in each sample
  mutate(fraction = counts/sum(counts)) %>%
  #log2 fold change of each barcode in each sample before and after treatment
  mutate(log2FC = log2(counts/n0)) %>%
  #strain varience
  mutate(strain_variance = (1/(1+counts)+1/(1+n0))/log(4, base = exp(1)) ) %>%
  #the gene fitness is the weighted average of the strain fitness
  #calculate the weight of each barcode
  mutate(weight_cal = 1/strain_variance) %>%
  #set a ceiling on the weight, being that of a strain with 20 reads in each sample
  mutate(weight_max = 1/(1/(1+20)+1/(1+n0))/log(4, base = exp(1))) %>%
  #the weight = min (weight_cal, weight_max)
  mutate(weight = ifelse(weight_cal < weight_max, weight_cal, weight_max)) %>%
  group_by(locusId) %>%
  #the number of strains in each gene (the number of barcodes inserted in each gene)
  mutate(strains_per_gene = n_distinct(barcode)) %>%
  ungroup() %>%
  #the sum of counts and n0 over all strains in all genes
  mutate(sum_nafter = sum(counts)) %>%
  mutate(sum_n0 = sum(n0))

#calculate gene fitness for genes with different strains separately
#for genes with only 1 or 2 strains:
genes_1or2_strains_k = filter(gene_fitness_k, strains_per_gene < 3) %>%
  mutate(pseudocounts = sum_nafter/sum_n0) %>%
  mutate(strain_fitness = log2(counts + sqrt(pseudocounts)) - log2(n0 + sqrt(1/pseudocounts))) %>%
  group_by(sample, locusId) %>%
  mutate(gene_fitness = sum(weight * strain_fitness)/sum(weight))

#for genes with more than 2 strains:
genes_3_strains_k = filter(gene_fitness_k, strains_per_gene >= 3) %>%
  group_by(sample, locusId) %>%
  mutate(pre_strain_fitness = log2(counts + 1) - log2(n0 + 1)) %>%
  mutate(pre_gene_fitness = median(pre_strain_fitness)) %>%
  mutate(pseudocounts = (2^pre_gene_fitness) * (sum_nafter/sum_n0)) %>%
  mutate(strain_fitness = log2(counts + sqrt(pseudocounts)) - log2(n0 + sqrt(1/pseudocounts))) %>%
  mutate(gene_fitness = sum(weight * strain_fitness)/sum(weight))

#concatenate the two tables of genes with different number of strains:
gene_fitness_k = bind_rows(select(genes_1or2_strains_k, -weight_cal, -weight_max, -sum_nafter, -sum_n0, -pseudocounts),
                           select(genes_3_strains_k, -weight_cal, -weight_max, -sum_nafter, -sum_n0, -pseudocounts, -pre_strain_fitness, -pre_gene_fitness))%>%
  mutate(identifier = "kyle")

# C. Correlation between qi's data & kyle's data:
#bind together the strain fitness scores of the same barcode from the same sample in kyle's and qi's data
barcode_bind = bind_rows(gene_fitness_qi, gene_fitness_k)

# (a) Correlation of strain fitness between kyle's data and qi's data for each sample.
cor_persample = barcode_bind %>%
  ungroup() %>%
  select(barcode, sample, strain_fitness, identifier) %>%
  #remove dupicate rows (some genes locate inside another gene, so some barcodes are inserted into two different genes)
  distinct() %>%
  group_by(identifier, sample, barcode) %>%
  mutate(strain_fitness_mean = mean(strain_fitness)) %>%
  select(-strain_fitness) %>%
  distinct() %>%
  spread(identifier, strain_fitness_mean) %>%
  #remove missing values and infinite values
  #then calculate the correlation of strain fitness between kyle's data and qi's data for each sample
  filter(is.finite(kyle) == TRUE & is.finite(qi) == TRUE) %>%
  group_by(sample) %>%
  summarise(correlation = cor(x = kyle, y = qi))
#save file
write_tsv(cor_persample, "correlation_per_sample")

# (b) Correlation between two replicates for qi's data and kyle's data respectively.
cor_replicates = barcode_bind %>%
  ungroup() %>%
  mutate(replicate = ifelse(endsWith(sample, "a"), yes = "A", no = "B")) %>%
  select(barcode, condition, strain_fitness, identifier, replicate) %>%
  #remove dupicate rows (some genes locate inside another gene, so some barcodes are inserted into two different genes)
  distinct() %>%
  group_by(identifier, condition, replicate, barcode) %>%
  mutate(strain_fitness_mean = mean(strain_fitness)) %>%
  select(-strain_fitness) %>%
  distinct() %>%
  spread(replicate,strain_fitness_mean) %>%
  filter(is.finite(A)==TRUE & is.finite(B)==TRUE) %>%
  group_by(condition, identifier) %>%
  summarise(correlation = cor(x = A, y = B)) %>%
  spread(identifier, correlation)
write_tsv(cor_replicates, "correlation_replicates")

# D. Principal component analysis:
# (a) PCA of qi's data (per sample):
#spread strain_fitness of each sample into columns
pca_qi = gene_fitness_qi %>%
  ungroup() %>%
  select(barcode, sample, strain_fitness) %>%
  distinct() %>%
  group_by(barcode, sample) %>%
  summarise(strain_fitness = mean(strain_fitness)) %>%
  spread(sample, strain_fitness) %>%
  ungroup()

#create a matrix
pcam_qi =  select(pca_qi, -barcode) %>%
  as.matrix()
rownames(pcam_qi) = pca_qi$barcode

#remove NA and Inf.
pcam_qi = pcam_qi[is.finite(rowSums(pcam_qi)), ] %>%
  scale(scale = FALSE)

PCA_qi = prcomp(pcam_qi)

#create plotting dataframe
PCAd_qi = as.data.frame(PCA_qi$rotation)
PCAd_qi$sample = rownames(PCAd_qi)
PCAd_qi = PCAd_qi %>%
  as_tibble() %>%
  #add date
  inner_join(select(barcode_qi %>% ungroup(), date, sample) %>% distinct())

library(scales)
library(ggrepel)

gp_qi = ggplot(PCAd_qi, aes(x=PC1, y=PC2, colour=date, label=sample))+
  geom_point(aes(size=PC3))+
  scale_size(range=c(1,3))+
  scale_colour_viridis_d()+
  geom_text_repel(force=3, size=3)
var_pc = percent(PCA_qi$sdev^2 / sum(PCA_qi$sdev^2))[1:3]
gp_qi = gp_qi+
  labs(
    x=paste("PC1 (", var_pc[1], ")", sep=""),
    y=paste("PC2 (", var_pc[2], ")", sep=""),
    size=paste("PC3 (", var_pc[3], ")", sep="")
  )+
  theme_bw()+
  theme(aspect.ratio=1)
gp_qi

ggsave(
  "pca_qi_strain_fitness.pdf", gp_qi,
  height=15/2.54, width=18/2.54
)

# (b) Strain_fitness distribution in each sample of qi' data:
pcat_qi = pcam_qi %>%
  as_tibble() %>%
  mutate(barcode = rownames(pcam_qi)) %>%
  gather(sample, strain_fitness, -barcode)

fs_qi = ggplot(pcat_qi, aes(x = strain_fitness)) +
  geom_histogram()+
  facet_wrap(~sample)
fs_qi

ggsave(
  "fs_distribution_qi.pdf", fs_qi,
  height=15/2.54, width=18/2.54
)

# (c) PCA of qi's data (per condition):
#in order to have more clear clusters:
pca_qi_2 = gene_fitness_qi %>%
  ungroup() %>%
  select(barcode, condition, sample, strain_fitness) %>%
  distinct() %>%
  #for each condition, use the mean strain fitness values of the two replicates
  group_by(barcode, condition) %>%
  mutate(strain_fitness_mean = mean(strain_fitness)) %>%
  select(-sample, -strain_fitness) %>%
  distinct() %>%
  spread(condition, strain_fitness_mean) %>%
  #remove abnormal samples 3a & 3b (Format)
  #remove time-zero samples (t0a, t0b, t02a, t02b, t03a, t03b)
  select(-Formate, -"Time0_8/14", -"Time0_8/16", -"Time0_8/28") %>%
  ungroup()

#create a matrix
pcam_qi_2 =  select(pca_qi_2, -barcode) %>%
  as.matrix()
rownames(pcam_qi_2) = pca_qi_2$barcode

#remove NA and Inf.
pcam_qi_2 = pcam_qi_2[is.finite(rowSums(pcam_qi_2)), ] %>%
  scale(scale = FALSE)

PCA_qi_2 = prcomp(pcam_qi_2)

#create plotting dataframe
PCAd_qi_2 = as.data.frame(PCA_qi_2$rotation)
PCAd_qi_2$condition = rownames(PCAd_qi_2)
PCAd_qi_2 = PCAd_qi_2 %>%
  as_tibble() %>%
  #add date
  inner_join(select(barcode_qi %>% ungroup(), date, condition) %>% distinct())

gp_qi_2 = ggplot(PCAd_qi_2, aes(x=PC1, y=PC2, colour=date, label=condition))+
  geom_point(aes(size=PC3))+
  scale_size(range=c(1,3))+
  scale_colour_viridis_d()+
  geom_text_repel(force=3, size=3)
var_pc = percent(PCA_qi_2$sdev^2 / sum(PCA_qi_2$sdev^2))[1:3]
gp_qi_2 = gp_qi_2+
  labs(
    x=paste("PC1 (", var_pc[1], ")", sep=""),
    y=paste("PC2 (", var_pc[2], ")", sep=""),
    size=paste("PC3 (", var_pc[3], ")", sep="")
  )+
  theme_bw()+
  theme(aspect.ratio=1)
gp_qi_2

ggsave(
  "pca_qi_2_strain_fitness.pdf", gp_qi_2,
  height=15/2.54, width=18/2.54
)

# (d) PCA of kyle's data:
#spread strain_fitness of each sample into columns
pca_k = gene_fitness_k %>%
  ungroup() %>%
  select(barcode, sample, strain_fitness) %>%
  distinct() %>%
  group_by(barcode, sample) %>%
  summarise(strain_fitness = mean(strain_fitness)) %>%
  spread(sample, strain_fitness) %>%
  ungroup()

#create a matrix
pcam_k =  select(pca_k, -barcode) %>%
  as.matrix()
rownames(pcam_k) = pca_k$barcode

#remove NA and Inf.
pcam_k = pcam_k[is.finite(rowSums(pcam_k)), ]%>%
  scale(scale = FALSE)

PCA_k = prcomp(pcam_k)

#create plotting dataframe
PCAd_k = as.data.frame(PCA_k$rotation)
PCAd_k$sample = rownames(PCAd_k)
PCAd_k = PCAd_k%>%
  as_tibble()%>%
  #add date
  inner_join(select(barcode_kyle%>% ungroup(), date, sample)%>% distinct())

gp_k = ggplot(PCAd_k, aes(x=PC1, y=PC2, colour=date, label=sample))+
  geom_point(aes(size=PC3))+
  scale_size(range=c(1,3))+
  scale_colour_viridis_d()+
  geom_text_repel(force=3, size=3)
var_pc = percent(PCA_k$sdev^2 / sum(PCA_k$sdev^2))[1:3]
gp_k = gp_k+
  labs(
    x=paste("PC1 (", var_pc[1], ")", sep=""),
    y=paste("PC2 (", var_pc[2], ")", sep=""),
    size=paste("PC3 (", var_pc[3], ")", sep="")
  )+
  theme_bw()+
  theme(aspect.ratio=1)
gp_k

ggsave(
  "pca_k_strain_fitness.pdf", gp_k,
  height=15/2.54, width=18/2.54
)

# (e) Strain_fitness distribution in each sample of kyle' data:
pcat_k = pcam_k %>%
  as_tibble() %>%
  mutate(barcode = rownames(pcam_k)) %>%
  gather(sample, strain_fitness, -barcode)

fs_k = ggplot(pcat_k, aes(x = strain_fitness)) +
  geom_histogram()+
  facet_wrap(~sample)
fs_k

ggsave(
  "fs_distribution_k.pdf", fs_k,
  height=15/2.54, width=18/2.54
)

# E. Hierarchical clustering:
hc_qi = t_test_qi %>%
  select(locusId, sample, norm_fg) %>%
  distinct() %>%
  spread(sample, norm_fg) %>%
  ungroup()

#creat a matrix
hcm_qi = select(hc_qi, -locusId) %>%
  as.matrix()
rownames(hcm_qi) = hc_qi$locusId

#remove NA and Inf.
hcm_qi = hcm_qi[is.finite(rowSums(hcm_qi)), ] %>%
  scale(scale = FALSE)

#calculate distance
dist_qi = dist(hcm_qi, method = "euclidean")
#hclust
hcluster_qi = hclust(dist_qi, method = "ward.D")

#plot colored dendrogram
library(dendextend)
plot(color_branches(hcluster_qi, k = 12))

#cutree into 12 clusters
cluster_qi = cutree(hcluster_qi, k = 12)
#make a tibble with genes, clusters and their fitness scores in each sample
genes_cluster_qi = cluster_qi %>%
  as_tibble() %>%
  mutate(locusId = names(cluster_qi)) %>%
  mutate(cluster = value) %>%
  select(-value)
genes_cluster_qi = inner_join(genes_cluster_qi, hc_qi) %>%
  #remove time-zero samples
  select(-t0a, -t0b, -t02a, -t02b, -t03a, -t03b) %>%
  gather(sample, gene_fitness, -locusId, -cluster)
#add conditions; calculate the mean gene_fitness of two replicates in each condition
cluster_condition_qi = left_join(genes_cluster_qi,
                                 gene_fitness_qi %>% ungroup() %>% select(sample, condition) %>% distinct()) %>%
  select(-sample) %>%
  group_by(locusId, condition) %>%
  mutate(gene_fitness_mean = mean(gene_fitness)) %>%
  select(-gene_fitness) %>%
  distinct()

#plot
# gp_cluster = ggplot(cluster_condition_qi, aes(x = condition, y = gene_fitness_mean, group = locusId)) +
#   xlab("Growth Conditions") + ylab("Gene Fitness") +
#   geom_line(aes(color = cluster)) +
#   geom_point(aes(color = cluster)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
#   theme(legend.position = "top") +
#   facet_wrap(~cluster, ncol = 2)
# gp_cluster

gp_cluster = ggplot(cluster_condition_qi, aes(x = condition, y = gene_fitness_mean, group = locusId)) +
  xlab("Growth Conditions") + ylab("Gene Fitness") +
  geom_line(aes(color = cluster)) +
  geom_point(aes(color = cluster)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.position = "top") +
  facet_grid(cluster~.)

ggsave(
  "clusters_condition.pdf", gp_cluster,
  height=(32/2.54), width=(18/2.54)
)

#plot of each cluster
gp_cluster_11 = ggplot(cluster_condition_qi %>% filter(cluster == "11"), aes(x = condition, y = gene_fitness_mean, group = locusId)) +
  xlab("Growth Conditions") + ylab("Gene Fitness") +
  geom_line() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  theme(legend.position = "top")
gp_cluster_11

# F. Gene Annotation:

Ralstonia = read_tsv("Ralstonia_annotations_by_Michael.tab")
genes_COG = cluster_condition_qi %>%
  ungroup() %>%
  select(locusId, cluster) %>%
  distinct() %>%
  inner_join(Ralstonia %>% select(COG_Process, locusId)) %>%
  group_by(cluster) %>%
  #the number of genes in each cluster
  mutate(sum = length(locusId)) %>%
  group_by(cluster, COG_Process) %>%
  #the number of each COG process in each cluster
  mutate(number = length(locusId)) %>%
  #the percentage of each COG process among all the processes in each cluster
  mutate(percent = number/sum) %>%
  select(-locusId, -sum) %>%
  distinct() %>%
  filter(COG_Process != "Function unknown") %>%
  #filter major COG processes (>= 10%) in each cluster
  filter(percent >= 0.1)
