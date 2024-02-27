# **Meeting #5**

## **Meeting Minutes** 
- **Should we include Westernized in our question?**
- **At what point do we merge the data sets? If we merge after before alpha and beta analysis will that affect when we're comparing the two cohorts?**
  
- **AIM 2: QIIME 2 Processing Results**
  
-   **What we have done so far**
-   created a merged_table.qza and merged_rep-seqs.qza file
-   created a merged_rooted-tree.qza
  
-  **Next steps**
-  training the classifier for the colombia dataset (in progress) - used the primers in the paper
-  training the classifier for the texas dataset - planning to use the V1 to V3 primers in this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8544895/ (Forward: 5'-3': AGA GTT TGA TYM TGG CTC AG Reverse: 5'-3': AGA GTT TGA TYM TGG CTC AG)

-  **Question:** Should we merge the classifiers and apply it the merged_rep_seqs.qza file OR should we apply the classifiers separately to each rep_seqs.qza file and merge the two taxonomy.qza files we generate?
-  filtering out non-bacterial DNA from our table

- **AIM 3 (up to 3C): Alpha and Beta Diversity Analysis Results**

- **General documentation:**
- Should we indicate if we used a detached screen to run certain jobs?
  
