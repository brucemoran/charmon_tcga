# **Cathal Harmon TCGA Survival Analysis**

## **Data and Analysis Overview**
Data was downloaded using functions provided by the [TCGAbiolinks R package](https://doi.org/10.1093/nar/gkv1507) for projects `TCGA-BRCA`, `TCGA-COAD`, `TCGA-LUAD`, `TCGA-OV` and `TCGA-UCEC`. Clinical survival data was downloaded from supplementary data provided in [Liu et al 2018](https://dx.doi.org/10.1016%2Fj.cell.2018.02.052). Survival analysis used our [rpartSurvivalClassifier]("https://github.com/DonaghEgan/rpartSurvivalClassifier") method based on the `rpart` package (Therneau, T., Atkinson, B., & Ripley, B. (2013). Rpart: Recursive Partitioning. R Package Version 4.1-3, http://CRAN.R-project.org/package=rpart) which allows stratification of the patient cohort based on gene expression and survival information. We also ran a conventional median gene expression stratification method. Survival analysis used the `survival` package (Therneau T (2021). A Package for Survival Analysis in R. R package version 3.2-13, https://CRAN.R-project.org/package=survival). We used the `survminer` package (https://github.com/kassambara/survminer) for visualisation.
