# Comparative Phylogenetic Comparisons of Tyrannus OTU Morphometrics
## Author: Maggie P. MacPherson
### Copyright: This documentation is available under a [Creative Commons]() licence


## Modification History
### See ['ReadMe Tyrannus morphology'](https://github.com/mmacphe/Tyrannus_morphology/commits/main/ReadMe_Tyrannus_morphology.Rmd)


## Purpose
### Prepare a phylogeny composed of all currently recognized Tyrannus subspecies with which to conduct comparative phylogenetic principal components analyses and phylogenetic ANOVA to assess morphological patterns across migratory, partially migratory, and sedentary lineages. 


## Preliminary Steps (not shown)
### 1. Download the Harvey et al. 2021 suboscine UCE phylogeny [here](https://github.com/mgharvey/tyranni#tyranni)
### 2. Extract the Tyrannus genus. This is the [MacPherson_Tyrannus_base.tre](https://github.com/mmacphe/Tyrannus_morphology/blob/main/MacPherson_Tyrannus_base.tre) file that is used in Step 1, below.


## Steps:
### 1. Follow [the 1_Tyrannus_phylogeny.R script](https://github.com/mmacphe/Tyrannus_morphology/blob/main/1_Tyrannus_phylogeny.R) to add in additional Tyrannus subspecies not present in the Harvey et al. 2021 UCE phylogeny. 

This creates the [Tyrannus_phylogeny.tre](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_phylogeny.tre) file used in subsequent comparative phylogenetic analyses.

### 2. Follow [the 2_morphometric_data_processing.R script](https://github.com/mmacphe/Tyrannus_morphology/blob/main/2_morphometric_data_processing.R) to process measurements from voucher specimens into the formats needed for comparative phylogenetic analyses. 

This code requires the following files as inputs: [Tyrannus_voucher_table.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Tyrannus_voucher_table.csv) [Tyrannus_phylogeny.tre](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_phylogeny.tre) (used to match OTU names in phylogeny to their names in the voucher table), and [Tyrannus_subspecies_MigrationStrategies.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Tyrannus_subspecies_MigrationStrategies.csv) (a list of OTUs and their assignment as migratory, partially migratory or sedentary). 

This code creates the following files as outputs: [Tyrannus morphology data.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus%20morphology%20data.csv) (morphometrics for all individuals along with names that match up with the names on the phylogeny), [morphology_summary_table.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/morphology_summary_table.csv) (the mean ± standard deviation and sample size of each morphometric for each OTU), [Tyrannus_data.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_data.csv) (the mean morphometric measurement for each OTU with names that match the names in the phylogeny and movement behaviour), and [cv_summary_table.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/cv_summary_table.csv) (the coefficient of variation for each morphometric for each OTU with names that match the names in the phylogeny and movement behaviour).

### 3. Follow [the 3_phylogenetic_PCA.R script](https://github.com/mmacphe/Tyrannus_morphology/blob/main/3_phylogenetic_PCA.R) to conduct phylogenetic principal components analysis of bill and feather morphometrics, and build a figure showing the phylogeny, ancestral state reconstruction of movement ecology and Bill PPC2, and PPCA plots for bills and 2 PPCA plots for feathers (one with all individuals, and one zoomed in to see the trend excluding most of the individuals from long-tailed OTUs). 

This code requires the following file as inputs: [Tyrannus_data.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_data.csv), [Tyrannus_phylogeny.tre](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_phylogeny.tre), [Tyrannus morphology data.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus%20morphology%20data.csv), [cv_summary_table.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/cv_summary_table.csv), [Tyrannus morphology + PCA avg.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus%20morphology%20%2B%20PCA%20avg.csv) (made during this part of the pipeline and used for th ancestral state reconstruction for the figure), and [Tyrannus_subspecies_MigrationStrategies.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Tyrannus_subspecies_MigrationStrategies.csv). 

This code creates the following files as outputs: [Tyrannus morphology + PCA avg.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus%20morphology%20%2B%20PCA%20avg.csv), [cv_summary.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/cv_summary.csv), [Tyrannus_phylogenetic_PCA.png](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_phylogenetic_PCA.png).

### 4. Follow [the 4_phyloANOVA.R script](https://github.com/mmacphe/Tyrannus_morphology/blob/main/4_phyloANOVA.R) to conduct the appropriate phylogenetic ANOVAs for each part of the dataset: morphometrics while accounting for body size (using tarsus length), PPCA scores, and coefficients of variation. 

This code requires the following files as inputs: [Tyrannus_phylogeny.tre](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus_phylogeny.tre), [Tyrannus morphology + PCA avg.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/Tyrannus%20morphology%20%2B%20PCA%20avg.csv), and  [cv_summary.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/cv_summary.csv). 

This code creates the following code as outputs: [phylANOVA_output.txt](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/phylANOVA_output.txt), [phylANOVA_output_PCscores.txt](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/phylANOVA_output_PCscores.txt), [phylANOVA_cv_output.txt](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/phylANOVA_cv_output.txt), and [phylANOVA_tarsus-corrected_residuals.csv](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Output%20Files/phylANOVA_tarsus-corrected_residuals.csv) (used to build boxplots in [Boxplots_Figures.R](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Boxplots_Figures.R)).



## Also included:
### 1. An .R script to produce our boxplot figures [here](https://github.com/mmacphe/Tyrannus_morphology/blob/main/Boxplots_Figures.R).

### 2. An .R script to assess age and sex class influences on mean morphometric values [here](https://github.com/mmacphe/Tyrannus_morphology/blob/main/LinearModels_Demographic_Influence.R). 



## Other things to note before diving right in to using this code: 
### 1. The .R scripts are all set up so that new runs save outputs to the source folder. All files created by the author are found in the [Output Files folder](https://github.com/mmacphe/Tyrannus_morphology/tree/main/Output%20Files).

### 2. Code to make separate output files, conduct separate analyses and create separate figures for females and males are commented out in the .R script files. Follow the instructions in each .R script to exchange lines of code with the sex-specific lines of code for analyses looking at sexes separately.

# Enjoy!
## Sincerely, 
Maggie

© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About


