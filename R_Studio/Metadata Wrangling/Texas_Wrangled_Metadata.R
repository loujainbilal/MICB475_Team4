####Texas_Wrangling####

#Load-in packages:
library(tidyverse)
library(phyloseq)
library(ape)
library(vegan) #just for visualizing the rarefaction plot
library(ggplot2)

#Export texas metadata
meta_texasFP <- "texas_metadata.tsv"
meta_texas <- read_delim(file=meta_texasFP, delim="\t")

###Format metadata table###

#Step 1: move "sample-id" to the row names.
samp_df <- as.data.frame(meta_texas[,-1])
rownames(samp_df) <- meta_texas$'sample-id'
SAMP <- sample_data(samp_df)

#Step 2: filter out all columns that arent needed. 
texas_filtered <- subset(SAMP, select = -c(barcode, proximal, distal, v_region,
                                           height, weight, sequence_run, 
                                           ALT, HOMA, experiment_center, sample_type, Description, 
                                           hip, arm_circ, pulse2, labid, 
                                           homa_ir, homa_beta))

#Step 3: Add and modify columns:

#Create a new vector for "prediabetic threshold"
#Empty new vector
texas_prediabetes <- c()

#Extract Yes/No's
for (x in texas_filtered$HbA1C) {
  if (x > 5.6) {
    texas_prediabetes <- c(texas_prediabetes, "Yes")
  } else {
    texas_prediabetes <- c(texas_prediabetes, "No")
  }
}

#Print the new vector
print(texas_prediabetes)

#Add prediabetes column:
texas_filtered$prediabetic <- texas_prediabetes

#Add "texas" column to keep track of datasets
texas_label <- c("Texas")
texas_filtered$dataset <-texas_label

#Change column names: 
colnames(texas_filtered)[1] <- "hba1c"
colnames(texas_filtered)[2] <- "sex"
colnames(texas_filtered)[3] <- "age"
colnames(texas_filtered)[4] <- "bmi"
colnames(texas_filtered)[5] <- "waist_circ"
colnames(texas_filtered)[6] <- "cholesterol"
colnames(texas_filtered)[7] <- "triglycerides"
colnames(texas_filtered)[8] <- "hdl"
colnames(texas_filtered)[9] <- "ldl"
colnames(texas_filtered)[10] <- "crp"
colnames(texas_filtered)[11] <- "fast_glucose"

#Change column order:
new_column_order2 <- c("hba1c", "fast_glucose", "age", "sex", "bmi", 
                       "crp", "cholesterol", "hdl", "ldl", "triglycerides", 
                       "waist_circ", "prediabetic", "dataset")
texas_filtered_final <- texas_filtered[, new_column_order2]

#Additional column with both prediabetic and dataset:
texas_filtered_final$prediabetic_and_dataset <- paste(texas_filtered_final$prediabetic, 
                                                      texas_filtered_final$dataset, sep = "_")

#Convert M/F from sex column:
sex_replacement <- c("male" = "M", "female" = "F")

texas_filtered_final$sex <- ifelse(texas_filtered_final$sex %in% names(sex_replacement), 
                                    sex_replacement[texas_filtered_final$sex], 
                                   texas_filtered_final$sex)

#Export texas_wrangled_metadata.tsv:

row_names <- rownames(texas_filtered_final)
column_names <- colnames(texas_filtered_final)
tsv_text <- apply(texas_filtered_final, 1, function(row) paste(row, collapse = "\t"))
tsv_text <- c(paste(column_names, collapse = "\t"), tsv_text)
tsv_text <- c(paste(c("", row_names), tsv_text, sep = "\t"))
writeLines(tsv_text, "texas_wrangled_metadata.tsv")


########A.S########
