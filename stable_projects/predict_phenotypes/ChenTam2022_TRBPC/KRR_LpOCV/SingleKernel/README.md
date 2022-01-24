# Kernel Regression Leave-p-Out Cross-Validation workflow

Functions in this folder provides a wrapper function to perform single-kernel leave-p-sites-out cross-validation workflow used in the paper. Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times.

# Usage
Run `./CBIG_TRBPC_KRR_LpOCV_workflow.sh --help` for the usage