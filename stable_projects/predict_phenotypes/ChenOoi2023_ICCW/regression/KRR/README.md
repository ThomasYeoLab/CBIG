# Kernel Regression Leave-p-Out Cross-Validation workflow

Functions in this folder provides a wrapper function to perform single-kernel leave-p-sites-out cross-validation workflow for KRR. Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times.

## Usage
Run the scripts in the following order: 
1. `CBIG_ICCW_KRR_LpOCV_prepare_parameters.m`: Prepare the folds and settings required for KRR
2. `CBIG_ICCW_run_rs_krr.sh`: Run KRR for specified sample size
3. `CBIG_ICCW_compute_krr_weights_job.sh`: Extract regression weights for KRR for specified sample size

Please read the instructions included in the scripts for directions on how to run the script. For shell scripts, you can enter the argument `--help` to see the usage.