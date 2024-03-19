####################################################
####Prediabetes wrangling final - March 19, 2024####
####################################################

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

###Create a new vector for "hba1c_prediabetic"###
#Empty new vector
hba1c_vector <- c()
#Extract Yes/No's for hba1c
for (x in col_filtered$Hemoglobin_a1c) {
  if (x >= 5.6) {
    hba1c_vector <- c(hba1c_vector, "Yes")
  } else {
    hba1c_vector <- c(hba1c_vector, "No")
  }
}
#Print the new vector
print(hba1c_vector)
#Add hba1c_prediabetes column:
col_filtered$hba1c_prediabetic <- hba1c_vector

###Create a new vector for "glucose_prediabetic"###
#Empty new vector
glucose_vector <- c()
#Extract Yes/No's for glucose
for (x in col_filtered$fast_glucose) {
  if (x >= 6.1) {
    glucose_vector <- c(glucose_vector, "Yes")
  } else {
    glucose_vector <- c(glucose_vector, "No")
  }
}
#Print the new vector
print(glucose_vector)
#Add glucose_prediabetic column:
col_filtered$fast_glucose_prediabetic <- glucose_vector

###Create a new vector for "hba1c_glucose_prediabetic"###
#Empty new vector
hba1c_glucose_vector <- c()
#Extract Yes/No's for hba1c_glucose
for (x in 1:nrow(col_filtered)) {
  if (col_filtered$Hemoglobin_a1c[x] >= 5.6 & col_filtered$fast_glucose[x] >= 6.1) {
    hba1c_glucose_vector <- c(hba1c_glucose_vector, "Yes")
  } else {
    hba1c_glucose_vector <- c(hba1c_glucose_vector, "No")
  }
}
#Print the new vector
print(hba1c_glucose_vector)
#Add hba1c_prediabetes column:
col_filtered$hba1c_glucose_prediabetic <- hba1c_glucose_vector

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
                      "waist_circ", "hba1c_prediabetic", "fast_glucose_prediabetic", "hba1c_glucose_prediabetic")
col_filtered_final <- col_filtered[, new_column_order]

#Convert M/F from sex column:
sex_replacement <- c("male" = "M", "female" = "F")

col_filtered_final$sex <- ifelse(col_filtered_final$sex %in% names(sex_replacement), 
                                   sex_replacement[col_filtered_final$sex], 
                                 col_filtered_final$sex)

#Export prediabetes_wrangled_metadata.tsv:

row_names <- rownames(col_filtered_final)
column_names <- colnames(col_filtered_final)
tsv_text <- apply(col_filtered_final, 1, function(row) paste(row, collapse = "\t"))
tsv_text <- c(paste(column_names, collapse = "\t"), tsv_text)
tsv_text <- c(paste(c("", row_names), tsv_text, sep = "\t"))
writeLines(tsv_text, "prediabetes_wrangled_metadata.tsv")

######A.S.######