# GEMPS
Biologically insightful survival analysis based on relating variability in gene expression phenotypes to genomic variants

The following packages are required:
- survival: for the Cox proportional hazards model
- randomForest: for imputing missing data
- survcomp: for computing the C-Index

The project has three files that you need to add to your current directory:
- Survival_Analysis_GEMPS.R: the main function and you need to source it
- Score: it will be sourced and called from the main function
- Predict: it will be sourced and called from the main function

Data should be imputed first if there are any missing values.

The output will be a list of two entries:
- The C-Index
- The indices of the gene expression features selected by GEMPS
