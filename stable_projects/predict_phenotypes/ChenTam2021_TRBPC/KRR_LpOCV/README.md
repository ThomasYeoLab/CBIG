# Kernel Ridge Regression Leave-p-Out Cross validation (KRR_LpOCV)

This folders contains leave-p-out cross validation workflow for both multi-kernel ridge regression and single-kernel ridge regression.
* `MultiKernel` contains the functions to run the LpOCV workflow for multi-kernel ridge regression
* `SingleKernel` contains the functions to run the LpOCV workflow for single-kernel ridge regression
# Leave-p-Out Cross validation
Suppose all the subjects are from N sites, in each cross-validation fold we choose subjects from p sites as test set and the remaining subjects in the training set. This is repeated N-choose-p times.