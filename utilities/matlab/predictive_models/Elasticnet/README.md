# Elastic-net Regression

This function uses glmnet to perform a cross-validated elastic-net regression.There are 2 parameters that are optimised for elasticnet: lambda (the regularisation parameter), and alpha (the ratio of L1 to L2 regularisation). Alpha and lambda are optimized through a gridsearch by choosing the parameters leading to the best prediction in the innerfold cross-validation in a given set of lambda and alpha.

Once calling the main function `CBIG_run_Elasticnet_workflow.m`, the regression is run in 3 main steps.
### Step 1:
The covariates are regressed from the target variable.
### Step 2.
The optimal hyperparameters are found in an inner-loop cross-validation in the function `CBIG_Elasticnet_innerloop_cv_glmnet.m`.
### Step 3:
The training accuracies for each fold of the outer-loop cross-validation is found by `CBIG_Elasticnet_train_test_glmnet.m`.

### Citation information for glmnet

Please use the following citation when using the glmnet package.

```
Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N.
http://www.stanford.edu/~hastie/glmnet_matlab/
```

----

## Usage

To run the elastic-net regression, a struct `params` with the following fields is used as an input to `CBIG_run_Elasticnet_workflow.m`.

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

   - `params.split_name` (compulsory)
     A string of the name of the split in the training test fold.

   - `params.outdir` (compulsory)
     A string of the full-path output directory.
 
   - `params.outstem` (compulsory)
     A string appended to the filename to specify the output files (y after regression, accuracy files, ...). For example, if `outstem = 'AngAffect_Unadj'`, then the final optimal accuracy file will be names as `fullfile(outdir, 'results', 'optimal_acc', 'AngAffect_Unadj.mat')`, and the optimal model hyperparameters of k-th fold will be saved in `fullfile(outdir, 'params', ['dis_' k '_cv'], 'selected_parameters_AngAffect_Unadj.mat')`.

   - `params.glmnet_dir` (optional)
     A string, the full path of the directory to store the glmnet code. Default is `$CBIG_CODE_DIR/external_packages/matlab/non_default_packages/glmnet/glmnet_matlab`.

   - `params.num_inner_folds` (optional)
     A string or a scalar. The number of inner-loop folds within each training fold for hyperparameter selection. If params does not have a field called `num_inner_folds`, this script uses the same number of outer-loop folds as the number of inner-loop folds.

   - `params.lambda` (optional)
     A vector of regularization hyperparameters to search through. The best lambda for each training-test split will be selected from the vector. The default set of vectors is [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20].

   - `params.alpha` (optional)
     A vector of ratios of L1 to L2 regularization to search through. The best alpha for each training-test split will be selected from the vector. The default set of vectors is [0.01 0.1 0.4 0.7 0.9 0.99].

   - `params.metric` (optional)
     A string, specifying how the value of lambda will be optimised with. Can choose the following: `corr`, `COD`, `predictive_COD`, `MAE`, `MAE_norm`,`MSE`,`MSE_norm`. Metric is `predictive_COD` by default.

   - `params.norm` (optional)
     A logical, specifying if features should be normalized. If true, features will be normalized over subjects in each test/training fold. False by default.


### Saved results

When `CBIG_run_Elasticnet_workflow.m` is run, a folder with the regressed y variables, `y`, is created in the output directory, containing the original y values and the y values after the covariates have been regressed for each fold. Additionally, a `params` folder is generated, where the optimal parameters for each fold is saved under `selected_parameters_<outstem>.mat` and the final accuracy and predicted y values from each fold of the inner loop is saved under `acc_train_<outstem>.mat`.

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

- Release v0.22.1 (08/03/2022): Fix bug in regressing covariates in functions calling `CBIG_regress_X_from_y_test.m`.
- Release v0.17.3 (13/01/2021): Initial release of Elasticnet regression package

----
## Bugs and questions
Please contact Leon Ooi at leonooiqr@gmail.com.
