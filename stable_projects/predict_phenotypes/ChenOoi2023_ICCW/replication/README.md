This folder contains scripts about how to replicate results of this project. Please ensure that the appropriate paths to the ABCD data have been exported from the config file using `CBIG_MMP_tested_config.sh`. The current paths are appropriate for those in the CBIG laboratory. Please modify these paths for your own environment.

# DATA
Input data used for replication is stored in `$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/ChenOoi2023_ICCW`. The data is contained in the `input` folder which contains the subject IDs, subject scores, behaviors to predict and covariates. The resting-FC matrix is stored in `$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/ChenTam2022_TRBPC` 
and not duplicated in this project's folder because of its large size.

# REPLICATION OF ANALYSIS
The analysis can be replicated in 3 steps. We provide 4 wrapper scripts to do so. The 4 wrapper scripts need to be run in a sequential manner as each step depends on the last. 
In each script, please modify `output_dir` to a directory of your choice. `output_dir` should be the same for all 4 wrapper scripts.

## Run prediction models
* First run `CBIG_ICCW_KRR_replication_wrapper.sh` to generate the KRR results.
* Next, run `CBIG_ICCW_LRR_LASSO_RF_replication_wrapper.sh` to generate the results for LRR, LASSO and RF. Note that the the scripts for the other regressions have been modified just to run for one behavior to save time. To replicate all behaviors please modify the scripts.

## Get predictive features from models and compare the agreement of predictive features
This step can only be run after the predictive models have finished running.
* First, run `CBIG_ICCW_interpretation_replication_wrapper.sh` which extracts the original regression weights from KRR, conditional variable importance from RF and t-statistics.
* Next, run `CBIG_ICCW_get_icc_replication_wrapper.sh` to calculate Haufe-transformed regression weights from all other models and find the agreement of predictive features.

## Results structure
After running the prediction models, in the output folder, there will be a folder for each regression containing the regression results and chosen hyperparameters for each fold and behavior. 
Additionally, once you run the interpretation and icc wrappers, there will be an analysis folder with the collated prediction accuracies and ICC scores for the original regression weights, Haufe-transformed regression weights and t-statistics.
A visualization of this structure is as follows

```
|-|output_folder
|---|KRR
|-----|800 #sample size
|-------|regression results
|-----|400 #sample size
|-------|regression results
|---|LRR / LASSO 
|-----|800 #sample size
|-------|behav_1 # behavior number
|---|RF 
|-----|800 #sample size
|-------|behav_1 # behavior number
|---------|fold_1 # fold number
|-----------|saved RF model
|-----|400 #sample size
|-------|regression results
|---|analysis
|-----|collated prediction accuracies and ICCs
```
