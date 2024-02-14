#Meeting #3 Notes:#
*Agenda*
Research Question: Is there a difference in the microbiome of Latin populations in Colombia compared to the United States as a predictor of pre-diabetes?
Pre-diabetes "yes" status: 
- HbA1C: 5.7 
- Fasting glucose: 100 to 125 mg/dL (5.6 to 6.9 mmol/L)

Questions
  - What's the point of using paper 7 if the data set is so small when considering HbA1C and fasting glucose?
      - Answer: Can we remove it
  - Paper 6 has values for ALL participants, Paper 7 only has values for diabetes, so technically we aren't just looking at pre-diabetes since one data set is only diabetic  
  - How to convert from TST or Excel --> TSV (Colombia data set) for data wrangling in R studio
      - Answer: can do this on numbers on a mac 
  - How do you use GitHub lol
  - "Diabetes" vs "glycemic control" vs "pre-diabetes"
      - Answer: combine pre-diabetic and diabetic (5.7 and above and 100 mg/dL and above) 
  - What's a good sample size for MICB 475?
      - Answer: 11 is small for Texas data set but we'll see where we're at after rarefaction 

*Meeting Minutes*
Texas
- Number of prediabetics:11
- Number of non-prediabetic: 52

Colombia
- Number of prediabetics: 135
- Number of non-prediabetic: 306

Date: 14th Feb 2023




**Meeting #2 Notes:**
Date: 7th Feb 2023

- Sources 6 and 7 had extensive metadata 
- Texas and Mexico data set in Project 2 on the server 
  - New data sets would be positive control and compare to Colombia's data set and see how similar their microbiome are 
- Filter Colombia column with pre-diabetic and not diabetic 
- See how the other two datasets defined diabetes/pre-diabetes
  - If then statements 
  - Download metadata from the server
    - Check for errors 
- How to classify diabetes?
- Two approaches to combining
  - Columbia - paired sequences
  - Paper number 6 - single-end sequences
  - Mexico - paired sequences
      - Look at the paper for the variable region
      - Process three separately into a table 
      - Combine into a table before diversity matrix 
1. Download all metadata 
2. All have a common column → on R
- One disease state column 
- To keep where they are from: Texas, Mexico, Columbia
3. Confounding variables 

- Datasets on the server: datasets→ project 2→ diabetes → Mexico and texas (both have metadata and manifest file)
- Proposal
    - Need to have the qiime2 pipeline done already
    - Get to the denoising step 
    - Don’t need to get all the metadata together

- Aims:
1. Data wrangling for metadata (metadata info)
- Got through all three, and have one column that has the diabetic state: the same way of defining 
- Remove what we don't need 
- A column with the same name and same definitions 
-Need R to merge metadata
2. Qiime2 processing (sequence processing and merging)
- Process them separately but will need to merge them later 
- Generate table.qzv files for each data 

- Limitations for the project:
    - Different data sets collected data differently, and processed data differently 
    - Maybe different metadata categories 

