###Load packages and files
#install the tidyverse
install.packages("tidyverse")

#install the ape package
install.packages("ape")

#install the phyloseq package 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

#load packages
library(tidyverse)
library(ape)
library(phyloseq)
library(vegan)
library(ggplot2)
library(microbiome)
library(ggVennDiagram)
library(DESeq2)

#install merged metadata 
metaFP <- "prediabetes_wrangled_metadata.tsv"
meta <- read_delim(file=metaFP, delim = "\t")

#install merged OTU table
otuFP <- "feature-table.txt"
otu <- read_delim(file=otuFP, delim = "\t", skip=1)

#install merged taxonomy file 
taxFP <-"taxonomy.tsv"
tax <- read_delim(file=taxFP, delim = "\t")

#install merged phylogenetic tree file 
treeFP <-"tree.nwk"
tree <- read.tree(file=treeFP)

###Adjust file format to be read by the phyloseq object
##format OTU table 
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$'#OTU ID'

#convert otu table to phyloseq readable object 
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

##format metadata
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df) <- meta$`sample-id`

#convert metadata to phyloseq readable object 
SAMP <- sample_data(samp_df)

##format taxonomy 
tax_mat <- tax %>% select(-Confidence)%>% 
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
#save and make sample ids the rows
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`

#make taxa table 
TAX <- tax_table(tax_mat)

#Create the phyloseq object 
team4 <- phyloseq(OTU, SAMP, TAX, tree)

#transform data to remove zeros 
team4_plus1 <- transform_sample_counts(team4, function(x) x+1)

#### HBA1C ###

#convert phyloseq object to deseq object: HBA1C 
team4_phyloseq_deseq_hba1c <- phyloseq_to_deseq2(team4_plus1, ~`hba1c_prediabetic`)

#run deseq command: HBA1C
DESEQ_team4_hba1c <- DESeq(team4_phyloseq_deseq_hba1c)
res_hba1c <- results(DESEQ_team4_hba1c, tidy=TRUE, 
                           #make No the reference group
                           contrast = c("hba1c_prediabetic","Yes","No"))
View(res_hba1c)

#### Fasting Glucose ###

#convert phyloseq object to deseq object: fasting glucose 
team4_phyloseq_deseq_fast_glucose <- phyloseq_to_deseq2(team4_plus1, ~`fast_glucose_prediabetic`)

#run deseq command: fasting glucose
DESEQ_team4_fast_glucose <- DESeq(team4_phyloseq_deseq_fast_glucose)
res_fast_glucose <- results(DESEQ_team4_fast_glucose, tidy=TRUE, 
                     #make No the reference group
                     contrast = c("fast_glucose_prediabetic","Yes","No"))
View(res_fast_glucose)

#### Both ###

#convert phyloseq object to deseq object: both 
team4_phyloseq_deseq_both <- phyloseq_to_deseq2(team4_plus1, ~`hba1c_glucose_prediabetic`)

#run deseq command: both
DESEQ_team4_both <- DESeq(team4_phyloseq_deseq_both)
res_both <- results(DESEQ_team4_both, tidy=TRUE, 
                            #make No the reference group
                            contrast = c("hba1c_glucose_prediabetic","Yes","No"))
View(res_both)

##### Visualizing results #####
##### Volcano plot #####

###### HBA1C #####

## Volcano plot: effect size VS significance

vol_plot_hba1c <- res_hba1c %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_hba1c

ggsave(filename="vol_plot_hba1c.png",vol_plot_hba1c)

# To get table of results

sigASVs_hba1c <- res_hba1c %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_hba1c)

# Get only asv names

sigASVs_hba1c_vec <- sigASVs_hba1c %>%
  pull(ASV)

sigASVs_hba1c_vec

###### Fasting glucose #####

## Volcano plot: effect size VS significance

vol_plot_fast_glucose <- res_fast_glucose %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_fast_glucose

ggsave(filename="vol_plot_fast_glucose.png",vol_plot_fast_glucose)

# To get table of results

sigASVs_fast_glucose <- res_fast_glucose %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_fast_glucose)

# Get only asv names

sigASVs_fast_glucose_vec <- sigASVs_fast_glucose %>%
  pull(ASV)

sigASVs_fast_glucose_vec

###### Both #####

## Volcano plot: effect size VS significance

vol_plot_both <- res_both %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_both

ggsave(filename="vol_plot_both.png",vol_plot_both)

# To get table of results

sigASVs_both <- res_both %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_both)

# Get only asv names

sigASVs_both_vec <- sigASVs_both %>%
  pull(ASV)

sigASVs_both_vec



