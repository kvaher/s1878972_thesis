# Bugs, brains and babies - The role of gut microbiota in preterm brain development

Code for the analyses presented in Chapter 5 and 6 of the PhD thesis

Repository author: Kadi Vaher (kadi.vaher@ed.ac.uk)

## Directory structure
The project assumes the following file directory system:

![image](https://user-images.githubusercontent.com/92927232/205872695-40a8f0e0-f4ea-498b-b03c-7ebc8b18e75c.png)

The repository consists of the following folders:
- **scripts/** containing R markdown/script files for creating the demographic/clinical variable descriptive tables, performing statistical analyses and creating accompanying figures and tables
- **src/** containing some custom functions used for loading data and creating figures

Note that metadata organisation scripts are not included in this repository.

## Chapter 5.	Microbiota profiles and drivers in preterm neonates
- microbiome_demotable.Rmd: code for creating study group demographic/clinical data descriptive tables
- qc_decontam.Rmd: code for the quality check and decontamination steps for the 16S rRNA sequencing data
- 16S_beta_diversity.Rmd: code to calculate Bray-Curtis dissimilarity matrix and conduct statistical analysis for clinical drivers of the beta diversity
- alpha_diversity.Rmd: code to calculate alpha diversity metrics and conduct statistical analysis for clinical drivers of the alpha diversity; also includes calculation of clinical variable differences between male and female infants
- hclust.Rmd: code for hierarchical clustering plot
- 16S_relative_abundance.Rmd: code for Maaslin2 models for the clinical covariates

## Chapter 6.	Neonatal microbiome and brain dysmaturation in preterm infants
- microbiota_mri_demo.Rmd: code to subset the participant list to be included in the microbiota-brain analyses; creating study group demographic/clinical data descriptive table
- 16S_ordination.Rmd: code to derive principal coordinates of microbiota data
- explore_braindata.Rmd: code to analyse the effects of GA at scan and GA at birth on MRI features of interest
- microbiota_brain_stats.Rmd: code for statistical analysis associating gut microbiota features with MRI features (regression and Maaslin2 models)
- dGM_regionaldiff_microbiome_stats.Rmd: code for statistical analysis associating gut microbiota features with MRI features in the deep grey matter regions



