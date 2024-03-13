
_**March 13th 2023**_

***Agenda***

**Taxa Bar Plot:**
![image](https://github.com/loujainbilal/MICB475_Team4/assets/159331304/ad4195b4-cdf7-45ef-b2b4-ce0ea60b61fc)

![image](https://github.com/loujainbilal/MICB475_Team4/assets/159331304/cc3122e8-396e-4b64-b81b-7969b2c7fc7c)


**Core microbiome results:**

Texas vs. Colombia

<img width="272" alt="Screenshot 2024-03-13 at 12 08 04 PM" src="https://github.com/loujainbilal/MICB475_Team4/assets/159101482/8d97b75f-1e7b-4d88-a535-0a0b38295f0e">


Prediabetic vs. Healthy

<img width="327" alt="Screenshot 2024-03-13 at 12 08 33 PM" src="https://github.com/loujainbilal/MICB475_Team4/assets/159101482/0a7ac669-b8c6-41b2-a10c-8074b50de443">


**Differential Abundance results:**

Prediabetic vs. Healthy

<img width="325" alt="Screenshot 2024-03-13 at 12 09 33 PM" src="https://github.com/loujainbilal/MICB475_Team4/assets/159101482/6476bd60-9f03-4aed-acf0-78213b5ac19c">

<img width="492" alt="Screenshot 2024-03-13 at 12 10 02 PM" src="https://github.com/loujainbilal/MICB475_Team4/assets/159101482/4ea03c13-add8-48d3-accc-85f3890f78df">


Texas vs. Colombia 

<img width="331" alt="Screenshot 2024-03-13 at 12 10 32 PM" src="https://github.com/loujainbilal/MICB475_Team4/assets/159101482/794e51e6-0f83-44d8-a74e-5cc784390e6f">


***Minutes***
- Taxabar plot, issues with assigning taxonomy to the Texas dataset, only eukaryote or unclassified. Can Simran please check over Aim 2 code?
    - yes something is wrong, we need to redo the classifier step
    - Do Texas and Colombia data sets separately
    - Can we merge separately classified files after they have been classified individually?
        - idk lol need to ask Dr. Sun or look online for code
    - Back up plan if assigning taxonomy to Texas doesn't work, just move forward with only Colombia dataset
    - If we end up only moving forward with Colombia dataset, we can look at other definitions of diabetes ( eg. include glucose)
      - Pick the most significant one and move forward with the downstream analysis
      - in the paper include the alpha and beta of all three categories (HbA1C, glucose, and both) ]
      - If none are statistically significant, use the one defined by the literature  
      - include functional analysis 
- 

***To Do***
- Pravin will assign taxa to two datasets SEPARATELY
- Merge before alpha-beta analysis 
  - Ask Dr. Sun if we can merge taxonomy files AFTER classified 
- Aiden to start backup meta-wrangling for just Colombia dataset 






