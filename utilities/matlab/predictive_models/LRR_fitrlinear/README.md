# Linear Ridge Regression using fitrlinear (LRR_fitrlinear)
This folder contains functions to perform linear ridge regression using matlab's built-in function fitrlinear as optimizer

# Usage
`CBIG_LRR_fitrlinear_workflow_1measure.m`
* This function is an adaptation to the CBIG_LRR_workflow_1measure in CBIG_repo. Instead of the usual inverse of matrices to fit the LRR, this function uses the fitrlinear function for fitting which is considerably faster. So instead of calling CBIG_LRR_train_test and CBIG_LRR_innerloop_cv, an adapted CBIG_LRR_fitrlinear_train_test and CBIG_LRR_fitrlinear_innerloop_cv is used.