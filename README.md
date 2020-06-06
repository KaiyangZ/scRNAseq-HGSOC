# scRNAseq-HGSOC
This repository contains the codes used in the paper "Analysis of single-cell RNA-seq data from ovarian cancer samples before and after chemotherapy links stress-related transcriptional profile with chemotherapy resistance". 

In this study, we developed a novel clustering method accounting for patient-specific variability and technical confounders, and we identified 12 clusters characterized by 10 distince gene communities from fresh tissue samples taken before and after neoadjuvant chemotherapy from 11 HGSOC patients.

### 1. Modeling and clustering scRNA-seq data
We modeled the observed single cell expression profiles as a mixture of latent transcriptional profiles and nuisance expression profiles following a Poisson distribution: 
                 **Y** ~ Poi( (( **X** %\*% **D** + **Z** %\*% **C** )) %\*% **G** )

where: **Y** are total expression profiles: total expression; **X** are nuisance; **D** is the design for nuisances: weights matrix; **Z** are latent expression profiles; **C** are latent classes (binary); **G** are sample scaling factors (e.g. total RNA). 

**Under `clustering/`**
* The alogrithm can be found in folder `poi_decom_gain` 
* In `fit_multk.scala`, we fitted our model with up to 25 clusters using 10 restarts, and the number of clusters is selected to be 12 based on Bayesian information criterion (BIC). 
* In `fit_k12.scala`, we performed 200 restarts for the selected number of clusters to ascertain a globally good clustering. 
* In `cluster2geneCommunities.R`, detection of co-expressed gene communities and subsequent analysis was performed with the clusters and profiles computed with the highest likelihood.

### 2. Stress-related transcriptional profile in TCGA data
* In `survival.R`, we checked the association between the progression-free survival and the detected stress-related transcriptional profile 





