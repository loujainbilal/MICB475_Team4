# **Meeting #6**

## **Meeting Minutes** 

- **What we worked on:**
-   finished the QIIME 2 processing - ended up not training a classifer and instead using the full sequences from the SILVA reference database
-   alpha and beta diversity core metrics using QIIME 2

- **Rarefaction** What are we supposed to base the sampling depth on? Dataset or just prediabetes? 
  - **Curve** 
  ![Screenshot 2024-03-05 at 12 27 46 AM](https://github.com/loujainbilal/MICB475_Team4/assets/159094203/680c9709-4f2c-4f0c-ae0c-0b333425dbe1)
![Screenshot 2024-03-05 at 12 27 54 AM](https://github.com/loujainbilal/MICB475_Team4/assets/159094203/605a0830-3d34-4253-a786-f7dbbf2e9c3b)
  
  - **Table**
![Screenshot 2024-03-05 at 12 24 01 AM](https://github.com/loujainbilal/MICB475_Team4/assets/159094203/778d44ed-cd29-4eb2-9270-a5ac0bae15cf)
![Screenshot 2024-03-05 at 12 24 13 AM](https://github.com/loujainbilal/MICB475_Team4/assets/159094203/0b0c3b1b-d76d-4ddd-b3f6-ab3635a11306)

- **Core diversity metrics** generated using rarefaction depth of 13890
  - n = 4 for yes texas
  - **Faiths pd**
    ![Screenshot 2024-03-05 at 12 16 35 AM](https://github.com/loujainbilal/MICB475_Team4/assets/159094203/216419c8-ec50-4f54-98ab-ce6a90de019c)

  - **Evenness**
    ![Screenshot 2024-03-05 at 12 19 10 AM](https://github.com/loujainbilal/MICB475_Team4/assets/159094203/388562fc-db7a-4ae8-b236-28f10f9688c5)


- **Notes**
  - If you plan on comparing between Colombia and Texas, the sampling depth should take all four categories into consideration
  - Sampling depth of 13890 is not great bc only retain 35% of features which is very low, but there is no other option as yes_texas sample size is very low and already at n = 4 for sampling depth of 13890
  - ideally would have done about 20000 to account for no_texas but then we would have less than n=3 for yes_texas
  - There are no barcodes because we merged two data sets, so rarefaction curve is grouping bases on clinical status or clinical status + data set

- **Action Items for Next Week**
  - Do the differential abundance first - Ayesha
      - depending on how this goes then we can see if we want to remove texas data set completely to use a higher sampling depth
  - Core microbiome analysis - Maya 
