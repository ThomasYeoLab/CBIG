This folder contains scripts about how to replicate results of this project. Notice that all filenames and directories in these scripts only work for CBIG lab.

# DATA
Input data used for replication is stored in `$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC`
# REPLICATE ALL ANALYSIS
`./CBIG_TRBPC_replication_all_wrapper.sh [outdir]`
* `outdir` is the output directory for the replication results.
* This function will run all the analysis done in the paper by calling all other scripts in this folder that only replicate part of our analysis.
* To skip some part of the analysis, just comment out the lines in this function.
* To check if the replication is correct, compare with the results in `$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/results`

# REPLICATE PART OF THE ANALYSIS
## Replicate the prediction models in the main analysis
`./CBIG_TRBPC_regression_main_wrapper.sh [top_outdir]`
* This function will run the regression models in the main analysis of our paper (kernel regression)
* `top_outdir` is the top level output directory. Inside this top level output directory, there will be sub-folders for each model. The folders will be:
  * `rs/` for single-kernel regression using resting-state functional connectivity (RSFC)
  * `mid/` for single-kernel regression using monetary incentive delay task FC
  * `nback/` for single-kernel regression using N-back task FC 
  * `sst/` for single-kernel regression using stop-signal task FC
  * `meanFC/` for single-kernel regression using average FC of all brain states. 
  * `allFC/` for multi-kernel regression using FC from all brain states.
## Replicate the prediction models in the control analysis
`./CBIG_TRBPC_regression_control_wrapper.sh [top_outdir]`
* This function will run the regression models in the main analysis of our paper (kernel regression with age and sex as covariates, and linear ridge regression)
* `top_outdir` is the top level output directory. Inside this top level output directory, there will be sub-folders for each model. The folders will be:
  * <top_outdir>/KRR_regress_age_sex/<feature_set> 
    * feature set is the features used in the model. It can be: `rs/`,`mid/`,`nback/`,`sst/`,`allFC/`
  * <top_outdir>/LRR/<feature_set>
    * for the linear ridge regression results. Instances in <feature_set> are the same as above.

## Replicate the predictive-feature matrices
* This funtion computes the predictive-feature matrix for all the regression models used in this paper. Usage: `CBIG_PFM_wrapper.sh [pred_result_dir] [PFM_dir]`
* `pred_result_dir`: the top level directory of prediction results. The same as `top_outdir` that you input to `CBIG_TRBPC_regression_main_wrapper.sh`
* `PFM_dir`: directory where you want to same the predictive feature matrix

## Replicate the permutation test for kernel regression
* This function runs permutation test for all the regression models used in the main analysis. Usage:`./CBIG_TRBPC_KRR_perm_test_wrapper.sh [pred_result_dir] [perm_outdir]`
* `pred_result_dir`: the top level directory of prediction results. The same as `top_outdir` that you input to `CBIG_TRBPC_regression_main_wrapper.sh`
* `perm_outdir`: output directory for the permutation results

## Replicate the permutation test for predictive-feature matrices
* This function runs permutation test for predictive-feature matrices of the multi-kernel FC model. Usage:`./CBIG_TRBPC_PFM_perm_test_wrapper.sh [pred_result_dir] [perm_outdir]`
* `pred_result_dir`: the top level directory of prediction results. The same as `top_outdir` that you input to `CBIG_TRBPC_regression_main_wrapper.sh`
* `perm_outdir`: output directory for the permutation results