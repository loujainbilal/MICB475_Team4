### Load packages ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library (picante)

#### Load data ####
metafp <- "prediabetes_wrangled_metadata_REDO.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "export_2.0/table_export_3.0/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "export_2.0/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "export_2.0/rooted-tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)


#### Format OTU table ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format metadata ####
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 

tax_mat <- tax_mat[,-1]

rownames(tax_mat) <- tax$`Feature ID`

TAX <- tax_table(tax_mat)
class(TAX)


#### Create phyloseq object ####
team4 <- phyloseq(OTU, SAMP, TAX, phylotree)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
team4_filt <- subset_taxa(team4,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
team4_filt_nolow <- filter_taxa(team4_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
team4_final <- prune_samples(sample_sums(team4_filt_nolow)>100, team4_filt_nolow)


# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
rarecurve(t(as.data.frame(otu_table(team4_final))), cex=0.1)
team4_rare <- rarefy_even_depth(team4_final, rngseed = 1, sample.size = 26493)


##### Saving #####
save(team4_final, file="team4_final.RData")
save(team4_rare, file="team4_rare.RData")


### ALPHA DIVERSITY ###

## RICHNESS ##

plot_richness(team4_rare) 

plot_richness(team4_rare, measures = c("Observed", "Shannon","Chao1")) 

#hba1c_prediabetic#
gg_richness_hba1c <- plot_richness(team4_rare, x = "hba1c_prediabetic", measures = c("Observed", "Shannon","Chao1")) +
  xlab("HbA1c Prediabetic Status") + 
  geom_boxplot(aes())
gg_richness_hba1c


ggsave(filename = "plot_three_hba1c.png"
       , gg_richness_hba1c
       , height=4, width=6)

#kruskal-wallis test for hba1c#

kruskal_sh_hba1c <- kruskal.test(Shannon ~ `hba1c_prediabetic`, data = samp_dat_wdiv)
# p-value = 0.08343, ns 

kruskal_chao1_hba1c <- kruskal.test(Chao1 ~ `hba1c_prediabetic`, data = samp_dat_wdiv)
# p-value = 0.1942, ns 

kruskal_observed_hba1c <- kruskal.test(Observed ~ `hba1c_prediabetic`, data = samp_dat_wdiv)
# p-value = 0.1659, ns 


#fast_glucose_prediabetic#
gg_richness_fpg <- plot_richness(team4_rare, x = "fast_glucose_prediabetic", measures = c("Observed", "Shannon","Chao1")) +
  xlab("Fasting Plasma Glucose (FPG) Prediabetic Status") +
  geom_boxplot()
gg_richness_fpg

ggsave(filename = "plot_three_fast_glucose.png"
       , gg_richness_fpg
       , height=4, width=6)


#fast_glucose_prediabetic#
gg_richness_fpg <- plot_richness(team4_rare, x = "fast_glucose_prediabetic", measures = c( "Shannon")) +
  xlab("Fasting Plasma Glucose (FPG) Prediabetic Status") +
  geom_boxplot()
gg_richness_fpg

ggsave(filename = "plot_shannon_fast_glucose.png"
       , gg_richness_fpg
       , height=4, width=3)

#kruskal-wallis test for fasting glucose#
kruskal_observed_fast <- kruskal.test(Observed ~ `fast_glucose_prediabetic`, data = samp_dat_wdiv)
# p-value = 0.5528, ns 

kruskal_chao1_fast <- kruskal.test(Chao1 ~ `fast_glucose_prediabetic`, data = samp_dat_wdiv)
# p-value = 0.6126, ns 

kruskal_sh_fast <- kruskal.test(Shannon ~ `fast_glucose_prediabetic`, data = samp_dat_wdiv)
# p-value = 0.35, ns 



## FAITH'S PD ##

#hba1c_prediabetic#
phylo_dist <- pd(t(otu_table(team4_rare)), phy_tree(team4_rare),
                 include.root=F) 

sample_data(team4_rare)$PD <- phylo_dist$PD

plot.pd <- ggplot(sample_data(team4_rare), aes(hba1c_prediabetic, PD)) + 
  geom_boxplot() +
  xlab("HbA1c Prediabetic Status") +
  ylab("Phylogenetic Diversity")

plot.pd

ggsave(filename = "pd_hba1c_prediabetic.png"
       , plot.pd
       , height=5, width=3)


#fast_glucose_prediabetic#
phylo_dist <- pd(t(otu_table(team4_rare)), phy_tree(team4_rare),
                 include.root=F) 

sample_data(team4_rare)$PD <- phylo_dist$PD

plot.pd <- ggplot(sample_data(team4_rare), aes(fast_glucose_prediabetic, PD)) + 
  geom_boxplot() +
  xlab("Fasting Glucose Prediabetic Status") +
  ylab("Phylogenetic Diversity")

plot.pd

ggsave(filename = "pd_glucose_prediabetic.png"
       , plot.pd
       , height=4, width=3)


#### BETA DIVERSITY - PERMANOVA ######

## JACCARD ##

dm_jaccard <- vegdist(t(otu_table(team4_rare)), method="jaccard")
adonis2(dm_jaccard ~ `hba1c_prediabetic`*`fast_glucose_prediabetic`, data=samp_dat_wdiv)

#                                             Df SumOfSqs      R2      F Pr(>F)
#hba1c_prediabetic                            1    0.422 0.00335 1.0352  0.338
#fast_glucose_prediabetic                     1    0.388 0.00308 0.9506  0.584
#hba1c_prediabetic:fast_glucose_prediabetic   1    0.369 0.00293 0.9051  0.835

## BRAY CURTIS ##

#hba1c_prediabetic#
bc_dm <- distance(team4_rare, method="bray")

pcoa_bc <- ordinate(team4_rare, method="PCoA", distance=bc_dm)

plot_ordination(team4_rare, pcoa_bc, color = "hba1c_prediabetic")

gg_pcoa_bc_hba1c <- plot_ordination(team4_rare, pcoa_bc, color = "hba1c_prediabetic") +
  labs(col = "HbA1c Prediabetic") + 
  stat_ellipse()
gg_pcoa_bc_hba1c

ggsave("plot_pcoa_bc_hba1c.png"
       , gg_pcoa_bc_hba1c
       , height=4, width=5)

#fast_glucose_prediabetic#
bc_dm <- distance(team4_rare, method="bray")

pcoa_bc <- ordinate(team4_rare, method="PCoA", distance=bc_dm)

plot_ordination(team4_rare, pcoa_bc, color = "fast_glucose_prediabetic")

gg_pcoa_bc_glucose <- plot_ordination(team4_rare, pcoa_bc, color = "fast_glucose_prediabetic") +
  labs(col = "Fasting Glucose Prediabetic") +
  stat_ellipse()
gg_pcoa_bc_glucose

ggsave("plot_pcoa_bc_glucose.png"
       , gg_pcoa_bc_glucose
       , height=4, width=5)


## Unweighted Unifrac ##

#hba1c_prediabetic#
unifrac_dis <- distance(team4_rare, method="unifrac")

pcoa <- ordinate(team4_rare, method="PCoA", distance = unifrac_dis)

plot_ordination(team4_rare, pcoa, color = "hba1c_prediabetic")

gg_pcoa_uf_hba1c <- plot_ordination(team4_rare, pcoa, color = "hba1c_prediabetic") +
  labs(col = "HbA1c Prediabetes") +
  stat_ellipse()
gg_pcoa_uf_hba1c 

## Saving plot ##
ggsave("plot_pcoa_uf_hba1c.png"
       , gg_pcoa_uf_hba1c
       , height=4, width=5)


#fast_glucose_prediabetic#
unifrac_dis <- distance(team4_rare, method="unifrac")

pcoa <- ordinate(team4_rare, method="PCoA", distance = unifrac_dis)

plot_ordination(team4_rare, pcoa, color = "fast_glucose_prediabetic")

gg_pcoa_uf_fast <- plot_ordination(team4_rare, pcoa, color = "fast_glucose_prediabetic") +
  labs(col = "Fasting Glucose Prediabetes") +
  stat_ellipse()
gg_pcoa_uf_fast 

## Saving plot ##
ggsave("plot_pcoa_uf_fast.png"
       , gg_pcoa_uf_fast
       , height=4, width=5)



## Weighted Unifrac ##

#hba1c_prediabetic#
wunifrac_dis <- distance(team4_rare, method="wunifrac")

pcoa <- ordinate(team4_rare, method="PCoA", distance = wunifrac_dis)

plot_ordination(team4_rare, pcoa, color = "hba1c_prediabetic")

gg_pcoa_wf_hba1c <- plot_ordination(team4_rare, pcoa, color = "hba1c_prediabetic") +
  labs(col = "HbA1c Prediabetes") +
  stat_ellipse()
gg_pcoa_wf_hba1c 

## Saving plot ##
ggsave("plot_pcoa_wf_hba1c.png"
       , gg_pcoa_wf_hba1c
       , height=4, width=5)

##WEIGHTED UNIFRAC - HbA1c##

dm_unifrac_hba1c <- UniFrac(team4_rare, weighted=TRUE)

ord.unifrac <- ordinate(team4_rare, method="PCoA", distance="unifrac")
plot_ordination(team4_rare, ord.unifrac, color="hba1c_prediabetic") +
  stat_ellipse(type = "norm")

adonis2(dm_unifrac_hba1c ~ `hba1c_prediabetic`*`fast_glucose_prediabetic`, data=samp_dat_wdiv)
#                                           Df SumOfSqs      R2      F Pr(>F)
#hba1c_prediabetic                            1 0.000795 0.00412 1.2695  0.256
#fast_glucose_prediabetic                     1 0.000267 0.00138 0.4267  0.820
#hba1c_prediabetic:fast_glucose_prediabetic   1 0.000358 0.00186 0.5719  0.495

## Saving plot ##
ggsave("plot_pcoa_wf_hba1c_2.png"
       , dm_unifrac_hba1c
       , height=4, width=5)

##WEIGHTED UNIFRAC - FPG##

dm_unifrac_fpg <- UniFrac(team4_rare, weighted=TRUE)

ord.unifrac <- ordinate(team4_rare, method="PCoA", distance="unifrac")
plot_ordination(team4_rare, ord.unifrac, color="fast_glucose_prediabetic") +
  stat_ellipse(type = "norm")

adonis2(dm_unifrac_fpg ~ `hba1c_prediabetic`*`fast_glucose_prediabetic`, data=samp_dat_wdiv)

#                                            Df SumOfSqs      R2      F Pr(>F)
#hba1c_prediabetic                            1 0.000795 0.00412 1.2695  0.288
#fast_glucose_prediabetic                     1 0.000267 0.00138 0.4267  0.778
#hba1c_prediabetic:fast_glucose_prediabetic   1 0.000358 0.00186 0.5719  0.526





#fast_glucose_prediabetic#
wunifrac_dis <- distance(team4_rare, method="wunifrac")

pcoa <- ordinate(team4_rare, method="PCoA", distance = unifrac_dis)

plot_ordination(team4_rare, pcoa, color = "fast_glucose_prediabetic")

gg_pcoa_wf_fast <- plot_ordination(team4_rare, pcoa, color = "fast_glucose_prediabetic") +
  labs(col = "Fasting Glucose Prediabetes") +
  stat_ellipse()
gg_pcoa_wf_fast 

## Saving plot ##
ggsave("plot_pcoa_wf_fast.png"
       , gg_pcoa_wf_fast
       , height=4, width=5)
