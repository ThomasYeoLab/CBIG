# Example of running prediction model and interpreting results.

This example shows how to run a kernel ridge regression (KRR) prediction model and extracts the predictive features by the Haufe transform.

## Data
Data for this example can be found in the `input` folder. This contains the data for 50 subjects. A description of the data for the subjects are as follows.
1. `y.mat`: A #subjects-long vector. The behavior measure to be predicted.
2. `RSFC.mat`: A 3D-matrix (#ROIs x #ROIs x #subjects). Contains the 40x40 symmetrical resting state-FC matrix
3. `covariates`: A #covariates x # subjects matrix. Contains 4 covariates to be regressed from the behavioral measures.
4. `no_relative_2_fold_sub_list.mat`: A struct containing 2 subfolds for the regression. Contains the subjects to be put in the training and test fold for the example.

## Scripts
1. `CBIG_ICCW_KRR_example_wrapper.m`: Runs the KRR regression and subsequently extracts the Haufe-transformed regression weights.
2. `CBIG_ICCW_check_example_results.m`: Checks whether the output of the example are the same as the reference results.