# LRR Leave-p-Out Cross-Validation workflow

Functions in this folder provides a wrapper function to perform single-kernel leave-p-sites-out cross-validation workflow for LRR. Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times. Note that alpha is set to 0 for Elasticnet, only L2 regularization is used, therefore it becomes linear ridge regression.

## Usage
Run the wrapper script `CBIG_ICCW_LRR_wrapper.sh`. Please ensure that KRR was run already, as these scripts requires the param files generated before running KRR. The wrapper script runs the following:
1. `CBIG_ICCW_LRR_get_parameters.m`: Copy param files from KRR and add additional required arguments for LRR. 
2. `CBIG_ICCW_Elasticnet_job.sh`: A shell script in the `regression/utilities` folder that calls the Elasticnet workflow.

Please read the instructions included in the scripts for directions on how to run the script. 