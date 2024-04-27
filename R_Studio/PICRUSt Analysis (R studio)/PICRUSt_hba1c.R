####################################################
####PICRUSt HbA1c - April 3, 2024####
####################################################

#Ensure library is loaded 
library(tidyverse)
library(ape)
library(phyloseq)
library(vegan)
library(ggplot2)

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

#Install PICRUST packages:
# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")
# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none
# After installing all of its above dependencies, install ggpicrust2
install.packages("ggpicrust2")
#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "pathway_abundance2.0.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  = as.data.frame(abundance_data)

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = meta$'sample-id'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names]

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
meta = meta[meta$`sample-id` %in% abun_samples,] #making sure the filtered metadata only includes these samples


####DESEQ ANALYSIS####


#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                        metadata = meta, group = "hba1c_prediabetic", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(421:ncol(abundance_desc))] 

# Generate a heatmap
pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = meta, group = "hba1c_prediabetic")

# Generate pathway PCA plot
pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = meta, group = "hba1c_prediabetic")

#Generate a Bar Graph
# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, meta, "hba1c_prediabetic")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05)
# You can also filter by Log2fold change

sig_res <- sig_res[order(sig_res$log2FoldChange),]
ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

####A.S###
