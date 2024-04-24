#Load-in packages:
library(tidyverse)
library(phyloseq)
library(ape)
library(vegan) #just for visualizing the rarefaction plot
library(ggplot2)


####Part 1: texas wrangling####


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


####Part 2: colombia wrangling####


#Export colombia metadata
meta_colFP <- "colombia_metadata.tsv"
meta_col <- read_delim(file=meta_colFP, delim="\t")

###Format metadata table###

#Step 1: move "#SampleID" to the row names.
samp2_df <- as.data.frame(meta_col[,-1])
rownames(samp2_df) <- meta_col$'#SampleID'
SAMP2 <- sample_data(samp2_df)

#Convert glucose units for colombia:
mg_dL_to_mmol_L <- function(mg_dL_value) {
  mmol_L_value <- mg_dL_value * 0.0555
  return(mmol_L_value)
}

SAMP2$fast_glucose <- mg_dL_to_mmol_L(SAMP2$glucose)

#Step 2: filter out all columns that arent needed. 
col_filtered <- subset(SAMP2, select = -c(adiponectin, age_range, BMI_class, Body_Fat_Percentage, 
                                          Calorie_intake, Cardiometabolic_status, city, country, diastolic_bp, 
                                          fiber, insulin, latitude, VLDL, medication, per_carbohydrates, 
                                          per_total_protein, per_total_fat, per_animal_protein, 
                                          per_monoinsaturated_fat, per_polyunsaturated_fat, per_saturated_fat, 
                                          smoker, stool_consistency, systolic_bp, MET_mins_per_week, glucose))

#Step 3: Add and modify columns:

#Create a new vector for "prediabetic threshold"
#Empty new vector
col_prediabetes <- c()

#Extract Yes/No's
for (x in col_filtered$Hemoglobin_a1c) {
  if (x > 5.6) {
    col_prediabetes <- c(col_prediabetes, "Yes")
  } else {
    col_prediabetes <- c(col_prediabetes, "No")
  }
}

#Print the new vector
print(col_prediabetes)

#Add prediabetes column:
col_filtered$prediabetic <- col_prediabetes

#Add "texas" column to keep track of datasets
col_label <- c("Colombia")
col_filtered$dataset <-col_label

#Change column names: 
print(colnames(col_filtered))

colnames(col_filtered)[1] <- "age"
colnames(col_filtered)[2] <- "bmi"
colnames(col_filtered)[3] <- "hba1c"
colnames(col_filtered)[4] <- "crp"
colnames(col_filtered)[5] <- "cholesterol"
colnames(col_filtered)[6] <- "hdl"
colnames(col_filtered)[7] <- "ldl"
colnames(col_filtered)[8] <- "triglycerides"
colnames(col_filtered)[9] <- "sex"
colnames(col_filtered)[10] <- "waist_circ"

#Change column order:
new_column_order <- c("hba1c", "fast_glucose", "age", "sex", "bmi", 
                      "crp", "cholesterol", "hdl", "ldl", "triglycerides", 
                      "waist_circ", "prediabetic", "dataset")
col_filtered_final <- col_filtered[, new_column_order]

 
####Part 3: merging and final formatting####


#Merge the two metadata files: 
team4_merged_metadata <- rbind(texas_filtered_final, col_filtered_final)

#Additional column with both prediabetic and dataset:
team4_merged_metadata$prediabetic_and_dataset <- paste(team4_merged_metadata$prediabetic, 
                                                       team4_merged_metadata$dataset, sep = "_")

#Convert M/F from sex column:
sex_replacement <- c("male" = "M", "female" = "F")

team4_merged_metadata$sex <- ifelse(team4_merged_metadata$sex %in% names(sex_replacement), 
                       sex_replacement[team4_merged_metadata$sex], 
                       team4_merged_metadata$sex)

#Export team4_metadata.tsv:

write.table(team4_merged_metadata, file = "team4_metadata.tsv", sep = "\t", row.names = TRUE, col.names = NA)

########A.S########
