####Colombia wrangling####

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

#Additional column with both prediabetic and dataset:
col_filtered_final$prediabetic_and_dataset <- paste(col_filtered_final$prediabetic, 
                                                    col_filtered_final$dataset, sep = "_")

#Convert M/F from sex column:
sex_replacement <- c("male" = "M", "female" = "F")

col_filtered_final$sex <- ifelse(col_filtered_final$sex %in% names(sex_replacement), 
                                   sex_replacement[col_filtered_final$sex], 
                                 col_filtered_final$sex)

#Export colombia_wrangled_metadata.tsv:

row_names <- rownames(col_filtered_final)
column_names <- colnames(col_filtered_final)
tsv_text <- apply(col_filtered_final, 1, function(row) paste(row, collapse = "\t"))
tsv_text <- c(paste(column_names, collapse = "\t"), tsv_text)
tsv_text <- c(paste(c("", row_names), tsv_text, sep = "\t"))
writeLines(tsv_text, "colombia_wrangled_metadata.tsv")



