# **Meeting #4**

## **Agenda** 
- AIMS!
   - Should we re-structure our aims because during our last meeting we established 6 aims which were general steps?
   - Ask Simran to look over "Data Wrangling" section outline to confirm contents.
        - F/U: were your suggestions in your email referring to the metadata data? Or the seq data? (we haven't yet looked into the seq data in great detail). 
        - F/U: We were going to only look at HBA1C or fasting glucose, not both, as many sources have noted that only looking at one is sufficient. We went with HBA1C as there is more literature. Does this work?
- Should we do our entire analysis using both the datasets from the start or do the analysis with one set (ex. Colombia), generate a model and apply it to the second set (ex. Texas)?
- Texas study found no significant associations with any of their clinical measures with any particular taxa? is this an issue for our study?
  - alpha beta diversity will determine this early on and just because they didnt find any significance doesn't mean we won't since we are combing with the Colombia data set so its larger 
- Colombia has manifest.txt file but Texas does not -> implications for importing step


## **Meeting Notes** ##

- AIMS
   - Restructure aims so they are similar to the proposal examples
- Aim 1: Data Wrangling
   -  does this need to be an aim? --> something vague and doesn't need to have an output that is being published, structure, challenges, and references for how we are defining pre-diabetes
   -  don't include specific limitations of variable regions - keep it to the metadata 
- Aim 2: QIIME2
   - QIIME2 goes into other aims - we'll have to process two different data separately before we can merge
   - End sequencing
      - Colombia is paired? but videos say all data sets are single end?
      - it shouldn't be paired - all course data sets are single end
   - Colombia has manifest.txt file but Texas does not -> implications for importing step
        - manifest tells you where all the folder are that you need, might need to create a manifest file for Texas, simran will try and look for one from published data and will check with Dr. Sun if not 
- Aim 3: Alpha and Beta Diversity
   - What groups are we comparing?
        - Colombia: prediabteic vs non diabetic
        - Texas: prediabteic vs non diabetic
    - pre vs not and whether col or USA, just compare pre vs not pre --> informs the downstream analysis 
        - should include alpha and beta diversity its one chunk of code 
        - might be affected by diet and lifestyle
        - if there is a difference in one or neither or both in either of the cohorts
        - okay to name taxa but not required  
- Aim 4: Indicator Taxa Analysis
   - no questions
- Aim 5: Core Microbiome Analysis
   - comparing pre vs not for now but if there is a significance in US vs Colombia then we could also do core microbiome for that  
- Aim 6: Differential Abundance
   - no questions
