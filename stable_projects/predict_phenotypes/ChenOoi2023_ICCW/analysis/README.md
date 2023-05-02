# Compare accuracy and ICC across models

Functions in this folder provides functions to: 
1. Collate the final accuracy (average over the two split-halves) 
2. Calculate the ICC of regression weights and Haufe transformed regression weights 
3. Calculate the ICC of the mass-univariate brain behavior relationships
4. Calculate the similarity between the regression weights and Haufe transformed regression weights 

# Usage
There is a folder for each regression model `KRR`, `LRR_LASSO`, `RF` in order to collate the accuracy and ICC for each of the models. Note that LRR and LASSO depend on the same script.
The `all_models` folder contains scripts to extract the mass-univariate brain behavior relationships (which is the same across models) and calculate the similarity of features extracted between models. 
For the latter, please ensure that the accuracy and ICC values have already been calculated. 

To collate the accuracy and ICC values for original regression weights and Haufe transformed weights for each model please run the following scripts:
1. `KRR/CBIG_ICCW_compute_KRR_acc_icc.m`: Collate results for KRR
2. `LRR_LASSO/CBIG_ICCW_compute_Elasticnet_acc_icc.m`: Collate results for LRR and LASSO. Please note that the LRR and LASSO depend on the same Elasticnet workflow as they are simply 2 cases of Elasticnet with alpha set to 0 and 1. 
3. `RF/CBIG_ICCW_compute_RF_acc_icc.m`: Collate results for RF

To calculate the mass-unvariate brain behavior relationships between FC and behaviors, use the following script:
1. `all_models/CBIG_ICCW_compute_tstats.m`: Collate mass univariate brain behavior relationship. Please note that the mass univariate brain behavior relationships are the same across all regression models, so we only have to calculate it once. In this case, we used the KRR results directory to do the calculation.

To calculate the similarity between the Haufe transformed / original regression weights, use the following two scripts:
1. `all_models/CBIG_ICCW_haufe_similarity.m`: Compare similarity for Haufe transformed weights across KRR, LRR, LASSO and RF
2. `all_models/CBIG_ICCW_weights_similarity.m`: Compare similarity for regression weights across KRR, LRR, LASSO (RF is not compared as it uses conditional variable importance)

Please read the instructions included in the scripts for directions on how to run the script. 