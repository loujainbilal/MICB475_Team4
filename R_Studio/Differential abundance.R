
#Differential Abundance analysis 

#load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)

#install merged metadata 
metaFP <- "team4_metadata_2.0.tsv"
meta <- read_delim(file=metaFP, delim = "\t")

#install merged OTU table
otuFP <- "feature-table.txt"
otu <- read_delim(file=otuFP, delim = "\t", skip=1)

#install merged taxonomy file 
taxFP <-"taxonomy.tsv"
tax <- read_delim(file=taxFP, delim = "\t")

#install merged phylogenetic tree file 
treeFP <-"tree.nwk"
tree <- read_tree(treeFP)

###Adjust file format to be read by the phyloseq object
##format OTU table 
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$'#OTU ID'

#convert otu table to phyloseq readable object 
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

##format metadata
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df) <- meta$`sample_name`

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

##Prediabetic comparison

#convert phyloseq object to deseq object: prediabetic
team4_phyloseq_deseq_prediabetic <- phyloseq_to_deseq2(team4_plus1, ~`prediabetic`)

#run deseq command: prediabetes
DESEQ_team4_prediabetes <- DESeq(team4_phyloseq_deseq_prediabetic)
res_prediabetes <- results(DESEQ_team4_prediabetes, tidy=TRUE, 
               #make No_columbia group the reference group
               contrast = c("prediabetic","Yes","No"))
View(res_prediabetes)

##Dataset comparison

#convert phyloseq object to deseq object: dataset
team4_phyloseq_deseq_dataset <- phyloseq_to_deseq2(team4_plus1, ~`dataset`)

#run deseq command: dataset
DESEQ_team4_dataset <- DESeq(team4_phyloseq_deseq_dataset)
res_dataset <- results(DESEQ_team4_dataset, tidy=TRUE, 
                           #make No_columbia group the reference group
                           contrast = c("dataset","Texas","Colombia"))
View(res_dataset)

### Visualizing results 

## Prediabetes 

## Volcano plot: effect size VS significance

vol_plot_prediabetes <- res_prediabetes %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_prediabetes

ggsave(filename="vol_plot_prediabetes.png",vol_plot_prediabetes)

# To get table of results

sigASVs_prediabetes <- res_prediabetes %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_prediabetes)

# Get only asv names

sigASVs_prediabetes_vec <- sigASVs_prediabetes %>%
  pull(ASV)

sigASVs_prediabetes_vec

## Dataset 

## Volcano plot: effect size VS significance

vol_plot_dataset <- res_dataset %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_dataset

ggsave(filename="vol_plot_dataset.png",vol_plot_dataset)

# To get table of results

sigASVs_dataset <- res_dataset %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_dataset)

# Get only asv names

sigASVs_dataset_vec <- sigASVs_dataset %>%
  pull(ASV)

sigASVs_dataset_vec


