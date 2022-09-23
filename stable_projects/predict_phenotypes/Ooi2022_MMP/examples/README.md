This folder contains a modified example of how to arrive at the results of these projects. The input data uses non-restricted data only. 

## DATA
Input data is stored in the `input` folder. The features provided are 400-length vector of cortical volume and the lower triange of the resting state FC matrix for 60 HCP subjects. For this example, we are performing a 3-fold cross validation, over a single split, without accounting for family structure. 

## Running the prediction
To run the example, run the first-level models first using `./CBIG_MMP_HCP_example_singlefeature_regression_wrapper.sh` first before running the second-level models. This is because the second-level models require the output from the first-level models to be used as features. The user needs to specify the output directory as an argument to the wrapper. 

### Single-feature prediction
To replicate the analysis, first run the first level models that predict behaviour using a single feature using the following wrapper:
`./CBIG_MMP_HCP_example_singlefeature_regression_wrapper.sh $output_dir`
* This function will predict the 3 factor scores in the HCP using the cortical volume and resting state features.
* KRR, LRR and Elasticnet regressions will be run, to skip any of these models, remove them from `$regression_arr` in the script.

## Multi-feature prediction
Secondly run the second level models that predict behaviour by combining multiple features. The output directory should be the same as the one from the first-level models.
`./CBIG_MMP_HCP_example_multifeature_regression_wrapper $output_dir`
* This function will predict behaviours for the 3 factor scores in the HCP. 
* Cortical volume and resting state features will be combined through multiKRR and stacking using LRR. To skip any of these models, remove them from `$regression_arr` in the script.

## Comparing to the reference
To check whether the output is correct, the user can call the following script to collate the results into a mat file.
`./examples/scripts/utilities/CBIG_MMP_HCP_collate_result_example_wrapper.m`
The script can then compare the difference to our tabulated results using the function `CBIG_MMP_check_example_results`. The user will need to specify the directory the results are kept, and whether he/she is checking the first or second level results.