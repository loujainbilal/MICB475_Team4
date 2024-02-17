# Meeting #3 Feb 14th Notes:
## *Agenda*

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

## *Meeting Minutes*
Texas
- Number of prediabetics:11
- Number of non-prediabetic: 52

Colombia
- Number of prediabetics: 135
- Number of non-prediabetic: 306

Date: 14th Feb 2023

General Notes: 
- Remove paper 7
- Focus on columbia and texas 
- Look at how people have defined pre-diabetic 
- Keep all columns when combining data sets and can filter out later 


Project Aims for proposal:

1. Data wrangling: To be done pre-proposal 
- Comprehensive data set with those two columns reconciled 
- Merge everything they have in common 
- Process them separately before making the phyloseq 
- Look into manifest + original paper to see which variable region they are looking et
- Will not have dataset overview section
- Call it meta data wrangling
- Explain what we did and how we classified pre-diabetic vs not
- Make third column to see where it came from (area)
  - Texas diabetic
  - Texas non-diabetic 
  - Columbia diabetic
  - Columbia non-diabetic 
- Dont remove rows

2. Qimme2 processing 
- Figure whether to combine from the beginning or after the table
  - Look at the regions where they sequences for 16S

3. Basic alpha and beta diversity 
- Comparing cohorts
- Comparing prediabetic to non diabetic 

4. Indicator taxa analysis 
5. Core microbiome analysis 
6. Differential abundance 
