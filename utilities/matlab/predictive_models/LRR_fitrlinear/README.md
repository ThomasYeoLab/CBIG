# Linear Ridge Regression using fitrlinear (LRR_fitrlinear)
This folder contains functions to perform linear ridge regression using matlab's built-in function fitrlinear as the optimizer.

# Usage
`CBIG_LRR_fitrlinear_workflow_1measure.m`
* This function is an adaptation to the CBIG_LRR_workflow_1measure in CBIG_repo. Instead of the usual inverse of matrices to fit the LRR, this function uses the fitrlinear function for fitting which is considerably faster. So instead of calling CBIG_LRR_train_test and CBIG_LRR_innerloop_cv, an adapted CBIG_LRR_fitrlinear_train_test and CBIG_LRR_fitrlinear_innerloop_cv is used.

To run the linear ridge regression, a struct `params` with the following fields is used as an input to `CBIG_LRR_fitrlinear_workflow_1measure.m`.

### Fields for `params`:
The struct should have the following fields:

   - `params.sub_fold` (compulsory)
     A structure containing the information of how the data are separated into training and test sets. `params.sub_fold(i).fold_index` is a #subjects x 1 binary vector. 1 refers to the corresponding subject is a test subject in the i-th test fold. 0 refers to the corresponding subject is a training subject in the i-th test fold.

   - `params.feature_mat` (compulsory)
     A matrix of the features used as the independent variable in the prediction. It can be a 2-D matrix with dimension of #features x #subjects or a 3-D matrix with dimension of #ROIs1 x #ROIs2 x #subjects. If params.feature_mat is 3-D, it is a connectivity matrix between two sets of ROIs, and only the lower-triangular off-diagonal entries will be considered as features because the connectivity matrix is symmetric.

   - `params.covariates` (compulsory)
     A matrix of the covariates to be regressed out from y (dimension: #subjects x #regressors). If the users only want to demean, an empty matrix should be passed in; if the users do not want to regress anything, set `params.covariates` = 'NONE'.

   - `params.y` (compulsory)
     A column vector of the target variable to be predicted (dimension: #subjects x 1).

   - `params.outdir` (compulsory)
     A string of the full-path output directory.
 
   - `params.outstem` (compulsory)
     A string appended to the filename to specify the output files (y after regression, accuracy files, ...). For example, if `outstem = 'AngAffect_Unadj'`, then the final optimal accuracy file will be names as `fullfile(outdir, 'results', 'optimal_acc', 'AngAffect_Unadj.mat')`, and the optimal model hyperparameters of k-th fold will be saved in `fullfile(outdir, 'params', ['dis_' k '_cv'], 'selected_parameters_AngAffect_Unadj.mat')`.

   - `params.num_inner_folds` (optional)
     A string or a scalar. The number of inner-loop folds within each training fold for hyperparameter selection. If params does not have a field called `num_inner_folds`, this script uses the same number of outer-loop folds as the number of inner-loop folds.

   - `params.gpso_dir` (optional)
     A string, the full path of the directory to store the cloned git repository for Gaussian process code. Current script needs the Gaussian-Process Surrogate Optimisation ([GPSO](https://github.com/jhadida/gpso)) package to optimize the objective function. This git repository and its dependency "deck" (https://github.com/jhadida/deck) need to be both cloned into the same folder, which is passed into current script through params.gpso_dir. Default is `$CBIG_CODE_DIR/external_packages/matlab/non_default_packages/Gaussian_Process`

   - `params.domain` (optional)
     The searching domain of parameters used by the Gaussian process optimization algorithm (dimension: #hyperparameters x 2). Default is [0.001 0.1; 3 8], where 0.001-0.1 is the searching domain for the feature selection threshold, and 3-8 is the searching domain for the L2 regularization hyperparameter (after taking logarithm). If feature selection is not needed a 1 x 2 vector with the search domain for lambda can be used. For more information, please refer to: https://github.com/jhadida/gpso

   - `params.eval` (optional)
     The maximal evaluation times of the objective function (a scalar). Default is 15. If it is too large, the runtime would be too long; if it is too small, the objective function may not be able to reach its optimal value. The users need to consider this trade-off when setting this variable. For more information, please refer to: https://github.com/jhadida/gpso
 
   - `params.tree` (optional)
     A scalar, the depth of the partition tree used to explore hyperparameters. Default is 3. For more information, please refer to: https://github.com/jhadida/gpso

   - `params.metric` (optional)
     A string, specifying how the value of lambda will be optimised with. Can choose the following: `corr`, `COD`, `predictive_COD`, `MAE`, `MAE_norm`,`MSE`,`MSE_norm`. Metric is `predictive_COD` by default.

### Saved results

When `CBIG_LRR_fitrlinear_workflow_1measure.m` is run, a folder with the regressed y variables, `y`, is created in the output directory, containing the original y values and the y values after the covariates have been regressed for each fold. Additionally, a `params` folder is generated, where the optimal parameters for each fold is saved under `selected_parameters_<outstem>.mat` and the final accuracy and predicted y values from each fold of the inner loop is saved under `acc_train_<outstem>.mat`.

The accuracies for the test folds are saved in `optimal_acc` under `<outstem>_final_acc.mat` with the following fields: 

   - `acc_corr_train`
     Training accuracy (given in correlation)

   - `acc_corr_test`
     Test accuracies (given in correlation)

   - `y_predict` 
     Predicted target values

   - `optimal_statistics`
     A cell array of size equal to the number of folds. Each element in the cell array is a structure storing the accuracies of each possible accuracy metric (eg. corr, MAE, etc).

----
## Updates

- Release v0.22.4 (31/03/2022): Add functionality to regress covariates from feature_mat.
- Release v0.22.1 (08/03/2022): Fix bug in regressing covariates in functions calling `CBIG_regress_X_from_y_test.m`.
- Release v0.17.3 (13/01/2021): Initial release of this fitrlinear-based linear ridge regression package.

----
## Bugs and questions
Please contact Jianzhong Chen at chenjianzhong1996@gmail.com or Leon Ooi at leonooiqr@gmail.com.
