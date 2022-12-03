# Bugs, brains and babies - The role of gut microbiota in preterm brain development

Code for the analyses presented in Chapter 5 and 6 of the PhD thesis

Repository author: Kadi Vaher (kadi.vaher@ed.ac.uk)

The repository consists of the following folders:
- **scripts/** containing R markdown and script files for creating the master dataframes, demographic and clinical variable tables, and statistical analyses 
- **src/** containing custom functions used for loading data and creating figures

## Chapter 5.	Microbiota profiles and drivers in preterm neonates
- qc_decontam.Rmd: code for the quality check and decontamination steps for the 16S rRNA sequencing data
- 16S_beta_diversity.Rmd: code to calculate Bray-Curtis dissimilarity matrix and conduct statistical analysis for clinical drivers of the beta diversity
- alpha_diversity.Rmd: code to calculate alpha diversity metrics and conduct statistical analysis for clinical drivers of the alpha diversity; also includes calculation of clinical variable differences between male and female infants
- hclust.Rmd: code for hierarchical clustering plot
- 16S_relative_abundance.Rmd: code for Maaslin2 models for the clinical covariates
- microbiome_demotable.Rmd: code for creating study group demographics descriptive tables
- sequencing_comparison.Rmd: code for the differences between meconium samples with and without sufficient bacterial yield for sequencing

## Chapter 6.	Neonatal microbiome and brain dysmaturation in preterm infants
- 16S_ordination.Rmd: code to derive principal coordinates of microbiota data
- explore_braindata.Rmd: code to analyse the effects of GA at scan and GA at birth on MRI features of interest
- microbiota_brain_stats.Rmd: code for statistical analysis associating gut microbiota features with MRI features (regression and Maaslin2 models)
- dGM_regionaldiff_microbiome_stats.Rmd: code for statistical analysis associating gut microbiota features with MRI features in the deep grey matter regions
