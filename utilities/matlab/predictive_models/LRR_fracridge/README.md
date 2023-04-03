# Linear Ridge Regression using fracridge (LRR_fracridge)
This folder contains functions to perform linear ridge regression using external toolbox https://github.com/nrdg/fracridge.

# Usage
`CBIG_LRR_frac_workflow.m`
* This function is an adaptation to the CBIG_LRR_workflow_1measure in CBIG_repo. Instead of the usual inverse of matrices to fit the LRR, this function uses the `fracridge` toolbox for optimization which is considerably faster. So instead of calling CBIG_LRR_train_test and CBIG_LRR_innerloop_cv, an adapted CBIG_LRR_frac_train_test and CBIG_LRR_frac_innerloop_cv is used.

Compared with other LRR scripts in our repo (i.e., LinearRidgeRegression and LRR_fitrlinear), one big advantage of this code is that the user can pass in multiple behavioral measures all at once instead of looping through each behavioral measure. Furthermore, the `fracridge` toolbox is extremely fast for data with large number of features (see Ariel Rokem, Kendrick Kay, Fractional ridge regression: a fast, interpretable reparameterization of ridge regression, GigaScience, Volume 9, Issue 12, December 2020, giaa133, https://doi.org/10.1093/gigascience/giaa133 for more details). In other LRR scripts in our repo, we will do feature selection (`CBIG_FC_FeatSel.m`) to speed up, which is sub-optimal. However, this is not necessary for this script because the `fracridge` toolbox is already very fast.

**IMPORTANT NODE:**
The current `fracridge` toolbox does not allow missing behavioral measures or missing features for participants. If there are missing behavioral measures or features for some participants, the script will check the features and target varibles and throw error messages if there are infinate values. 

To run the linear ridge regression, a struct `params` with the following fields is used as an input to `CBIG_LRR_frac_workflow.m`.

### Fields for `params`:
The struct should have the following fields:

   - `params.sub_fold` (compulsory)
     A structure containing the information of how the data are separated into training and test sets. `params.sub_fold(i).fold_index` is a #subjects x 1 binary vector. 1 refers to the corresponding subject is a test subject in the i-th test fold. 0 refers to the corresponding subject is a training subject in the i-th test fold.

   - `params.feature_mat` (compulsory)
     A matrix of the features used as the independent variable in the prediction. It can be a 2-D matrix with dimension of #features x #subjects or a 3-D matrix with dimension of #ROIs1 x #ROIs2 x #subjects. If params.feature_mat is 3-D, it is a connectivity matrix between two sets of ROIs, and only the lower-triangular off-diagonal entries will be considered as features because the connectivity matrix is symmetric.

   - `params.covariates` (compulsory)
     A matrix of the covariates to be regressed out from y (dimension: #subjects x #regressors). If the users only want to demean, an empty matrix should be passed in; if the users do not want to regress anything, set `params.covariates` = 'NONE'.

   - `params.y` (compulsory)
     A matrix of the target variables to be predicted (dimension: #subjects x K). K is the number of target variables to be predicted.

   - `params.lambda` (optional)
     A vector of the fraction parameters between 0 and 1. Fractions can be exactly 0 or exactly 1. However, values in between 0 and 1 should be no less than 0.001 and no greater than 0.999. For example, params.lambda could be 0:.05:1 or 0:.1:1.  If this field is not passed in, the script will use the default lambda values: params.lambda = [0.05:0.05:1]. 

   - `params.outdir` (compulsory)
     A string of the full-path output directory.
 
   - `params.outstem` (compulsory)
     A string appended to the filename to specify the output files (y after regression, accuracy files, ...). For example, if `outstem = 'AngAffect_Unadj'`, then the final optimal accuracy file will be named as `fullfile(outdir, 'results', 'optimal_acc', 'AngAffect_Unadj.mat')`, and the optimal model hyperparameters of k-th fold will be saved in `fullfile(outdir, 'params', ['dis_' k '_cv'], 'selected_parameters_AngAffect_Unadj.mat')`.

   - `params.num_inner_folds` (optional)
     A string or a scalar. The number of inner-loop folds within each training fold for hyperparameter selection. If params does not have a field called `num_inner_folds`, this script uses the same number of outer-loop folds as the number of inner-loop folds.

   - `params.metric` (optional)
     A string, specifying how the value of lambda will be optimised with. Can choose the following: `corr`, `COD`, `predictive_COD`, `MAE`, `MAE_norm`,`MSE`,`MSE_norm`. Metric is `predictive_COD` by default.

### Saved results

When `CBIG_LRR_frac_workflow.m` is run, a folder with the regressed y variables, `y`, is created in the output directory, containing the original y values and the y values after the covariates have been regressed for each fold. Additionally, a `params` folder is generated, where the optimal parameters for each fold is saved under `selected_parameters_<outstem>.mat` and the final accuracy and predicted y values from each fold of the inner loop is saved under `acc_train_<outstem>.mat`.

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

- Release v0.28.0 (03/04/2023): Initial release of LRR_fracridge. These scripts were used in Kong2023 RSFC prediction paper.

----
## Bugs and questions
Please contact Ruby Kong at roo.cone@gmail.com or Leon Ooi at leonooiqr@gmail.com.
