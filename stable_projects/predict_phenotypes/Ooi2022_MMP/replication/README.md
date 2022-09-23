This folder contains scripts about how to replicate results of this project. Please ensure that the appropriate paths to the ABCD and HCP data have been exported from the config file using `CBIG_MMP_tested_config.sh`. The current paths are appropriate for those in the CBIG laboratory. Please modify these paths for your own environment.

# DATA
Input data used for replication is stored in `$CBIG_HCP_REPDATA_DIR/data` and `$CBIG_ABCD_REPDATA_DIR/data` for the HCP and ABCD respectively. Each repository has a folder `behaviour` which contains the subject IDs, behaviour scores and covariates. Also, there will be a folder `features` which contain the imaging features derived from each modality. Each feature is saved in a #features x #subjects matrix. 

# REPLICATION OF ANALYSIS
We provide 8 scripts to replicate all aspects of this study (4 for the HCP and 4 for the ABCD). Replication scripts need to be run in a sequential manner as the stacking models depend on the output of the KRR models. **IMPORTANT: In each script, please modify `output_dir` to a directory of your choice.** 

## Single-feature prediction
To replicate the analysis, first run the first level models that predict behaviour using a single feature.
`./CBIG_MMP_HCP_singlefeature_regression_wrapper`
`./CBIG_MMP_ABCD_singlefeature_regression_wrapper`
* This function will predict behaviours for 58 behaviours in the HCP and 39 behaviours in the ABCD using each anatomical, diffusion and functional MRI feature presented in the paper.
* This script calls `Ooi2022_MMP/regression/<dataset>/utilities/CBIG_MMP_<dataset>_schedule_regression.sh` which contains the entire list of features from each modality. Modify the list from this script if you wish to run a subset of features.
* KRR, LRR and Elasticnet regressions will be run, to skip any of these models, remove them from `$regression_arr` in the script.
* To check if the replication is correct, compare with the results in `$CBIG_<dataset>_REPDATA_DIR/results`.

## Combined-feature prediction
Secondly run the second level models that predict behaviour by combining multiple features
`./CBIG_MMP_HCP_multifeature_regression_wrapper`
`./CBIG_MMP_ABCD_multifeature_regression_wrapper`
* This function will predict behaviours for 58 behaviours in the HCP and 39 behaviours in the ABCD using each anatomical, diffusion and functional MRI feature presented in the paper.
* This script calls `Ooi2022_MMP/regression/<dataset>/utilities/CBIG_MMP_<dataset>_schedule_regression.sh` which contains the entire list of features to combine. Modify the list from this script if you wish to change the features being combined.
* MultiKRR and stacking regressions using LRR will be run, to skip any of these models, remove them from `$regression_arr` in the script.
* To check if the replication is correct, compare with the results in `$CBIG_<dataset>_REPDATA_DIR/results`.

## Generate permutation models
Next, run the scripts that generate null models for the permutation test.
`./CBIG_MMP_HCP_perm_regression_wrapper`
`./CBIG_MMP_ABCD_perm_regression_wrapper`
* This function will permute the 3 factor scores among subjects for the HCP and ABCD to generate the null models.
* Null models will be created for all KRR, multiKRR and stacking models.

## Statistical tests
Lastly, the p-values and FDR correction will be performed for each dataset.
`./CBIG_MMP_HCP_perm_regression_wrapper`
`./CBIG_MMP_ABCD_perm_regression_wrapper`
* This function will permute the 3 factor scores among subjects for the HCP and ABCD to generate the null models.
* Null models will be created for all KRR, multiKRR and stacking models.