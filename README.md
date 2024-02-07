Meeting #2 Notes:

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

