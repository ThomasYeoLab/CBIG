# Linear Ridge Regression Leave-p-Out Cross validation (LRR_LpOCV)

This folders contains leave-p-out cross validation workflow for linear ridge regression using FC of single brain state or multiple brain-states.
# Leave-p-Out Cross validation
Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times.
# Usage
Run `./CBIG_TRBPC_LRR_LpOCV_workflow.sh --help` for the usage