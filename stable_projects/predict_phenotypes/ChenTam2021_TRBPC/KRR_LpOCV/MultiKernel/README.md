# Multi-Kernel Regression Leave-p-Out Cross-Validation workfolw

Functions in this folder provides a wrapper function to perform multi-kernel leave-p-sites-out cross-validation workflow used in the paper. Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times.

# Usage
`CBIG_TRBPC_multiKRR_LpOCV_workflow.sh` is the wrapper function that performs the  multi-kernel leave-p-sites-out cross-validation workflow.

Run `./CBIG_TRBPC_multiKRR_LpOCV_workflow.sh --help` for more information.