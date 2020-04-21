For mQTL analysis we need phenotype matrix (metabolites matrix) as normally distributed as possible.

Fist, we try to get rid of unwanted biases by using EigenMS script, then apply standard quantile normalization procedure.


All input files and R scripts to run this example are here: Files/Multimodal analyses/Normalised data 

Here is my normalization procedures for MS data:
# download both scripts into one folder (normalization_MS.R, EigenMS.R)
# from R console
source("normalization_MS.R")

1) 
EigenMS.  “Metabolomics Data Normalization with EigenMS”, Yuliya V. Karpievitch et al. 
EigenMS removed bias of unknown complexity from the LC-MS metabolomics data, allowing for increased sensitivity in differential analysis. EigenMS works in several stages. 
EigenMS preserves the treatment group differences in the metabolomics data by estimating treatment effects with an ANOVA model (multiple fixed effects can be estimated). 
Singular value decomposition of the residuals matrix is then used to determine bias trends in the data. 
The number of bias trends is then estimated via a permutation test and the effects of the bias trends are eliminated.

R function name: EigenMS_normalize(matrix_to_normalize,matrix_EigenMS_result,  covariates_to_preserve)

2) 
Quantile normalization
The transformation of the measurements for each metabolite into normally distributed while preserving relative rankings.

R function name:  QN_normalize(matrix_to_normalize,matrix_QN_result)


Example:
SHPOS_original.csv - original MS data
SHPOS_EigenMS.csv - file that will be created by R script
Covariates_EigenMS.txt - file with covariates that have to be preserved 
SHPOS_EigenMS_QN.csv - file that will be created by R script

1) EigenMS_normalize("..path../SHPOS_original.csv","..path_to../SHPOS_EigenMS.csv", "..path_to../Covariates_EigenMS.txt")

2) QN_normalize("..path_to../SHPOS_EigenMS.csv","..path_to../SHPOS_EigenMS_QN.csv")




