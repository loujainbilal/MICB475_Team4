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
metaFP <- "prediabetes_wrangled_metadata_REDO.tsv"
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
print(TAX)

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

###### HBA1C #####

## Volcano plot: effect size VS significance

#Color install:
install.packages("RColorBrewer")
library(RColorBrewer)
colors_blue <- brewer.pal(9, "Blues")

vol_plot_hba1c <- res_hba1c %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot_hba1c <- res_hba1c %>%
  mutate(significant = padj< 0.05 & abs(log2FoldChange)> 2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  scale_color_manual(values= c("#C6DBEF", "#08519C")) +
  labs(x = "Log2FoldChange", y = "-Log10(padj)") 

vol_plot_hba1c

ggsave(filename="vol_plot_hba1c_redo.png",vol_plot_hba1c)


# To get table of results

sigASVs_hba1c <- res_hba1c %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_hba1c)

write.csv(sigASVs_hba1c, file="Differential Abundance HBA1C.csv")

# Get only asv names

sigASVs_hba1c_vec <- sigASVs_hba1c %>%
  pull(ASV)

sigASVs_hba1c_vec

###### Fasting glucose #####

## Volcano plot: effect size VS significance

vol_plot_fast_glucose <- res_fast_glucose %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  scale_color_manual(values= c("#BCBDDC", "#54278F")) +
  labs(x = "Log2FoldChange", y = "-Log10(padj)") 


vol_plot_fast_glucose

ggsave(filename="vol_plot_fast_glucose_redo.png",vol_plot_fast_glucose)

# To get table of results

sigASVs_fast_glucose <- res_fast_glucose %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs_fast_glucose)

write.csv(sigASVs_fast_glucose, file="Differential Abundance Fasting Glucose.csv")

# Get only asv names

sigASVs_fast_glucose_vec <- sigASVs_fast_glucose %>%
  pull(ASV)

sigASVs_fast_glucose_vec

#### Prune phyloseq file: HBA1C ####

DESEQ_team4_hba1c <- prune_taxa(sigASVs_hba1c_vec,team4)

sigASVs_hba1c1 <- tax_table(DESEQ_team4_hba1c) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_hba1c) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs_hba1c1) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#### Prune phyloseq file: Fasting glucose ####

DESEQ_team4_fast_glucose <- prune_taxa(sigASVs_fast_glucose_vec,team4)

sigASVs_fast_glucose1 <- tax_table(DESEQ_team4_fast_glucose) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_fast_glucose) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs_fast_glucose1) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


##### INDICATOR TAXA ANALYSIS #####

#Install package
install.packages("indicspecies")

#Load packages
library(tidyverse)
library(phyloseq)
library(indicspecies)

# Remove non-bacterial sequences, if any
team4_filt <- subset_taxa(team4,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
team4_filt_nolow <- filter_taxa(team4_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
team4_final <- prune_samples(sample_sums(team4_filt_nolow)>100, team4_filt_nolow)

#Save RData File 
save(team4_final, file="team4_final_REDO.RData")

#Load Data
load("team4_final_REDO.RData")

# glom to Genus
team4_genus <- tax_glom(team4_final, "Genus", NArm = FALSE)
team4_genus_RA <- transform_sample_counts(team4_genus, fun=function(x) x/sum(x))

#ISA: HBA1C
isa_team4_hba1c <- multipatt(t(otu_table(team4_genus_RA)), cluster = sample_data(team4_genus_RA)$`hba1c_prediabetic`)
summary(isa_team4_hba1c)

taxtable_hba1c <- tax_table(team4_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")


# Create summary table 

isa_table_hba1c <- isa_team4_hba1c$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable_hba1c) %>%
  filter(p.value<0.05)

write.csv(isa_table_hba1c, file="ISA_HBA1C_REDO.csv")

#ISA: fasting glucose
isa_team4_fast_glucose <- multipatt(t(otu_table(team4_genus_RA)), cluster = sample_data(team4_genus_RA)$`fast_glucose_prediabetic`)
summary(isa_team4_fast_glucose)

taxtable_fast_glucose <- tax_table(team4_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Create summary table 

isa_table_fast_glucose <- isa_team4_fast_glucose$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable_fast_glucose) %>%
  filter(p.value<0.05)

write.csv(isa_table_fast_glucose, file="ISA_FASTING_GLUCOSE_REDO.csv")









