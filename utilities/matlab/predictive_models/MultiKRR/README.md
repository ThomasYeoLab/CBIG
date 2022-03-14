# Multi kernel ridge regression (MultiKRR)
This folder contains functions to perform the multi kernel ridge regression. 

# Usage
`CBIG_MultiKRR_workflow.m`
The multiKRR is an extension of the KRR workflow to perform behavior prediction while combining multiple feature types. To run the multiKRR regression, a struct `setup_file.mat` with the following fields is used as an input to `CBIG_MultiKRR_workflow.m`. Otherwise, the user can choose to input the arguments individually.

### Fields for `setup_file`:
The struct should have the following fields:

   - `sub_fold` (compulsory)
     A structure containing the information of how the data are separated into training and test sets. `sub_fold(i).fold_index` is a #subjects x 1 binary vector. 1 refers to the corresponding subject is a test subject in the i-th test fold. 0 refers to the corresponding subject is a training subject in the i-th test fold.

   - `num_inner_folds` (compulsory)
     A string or a scalar. The number of inner-loop folds within each training fold for hyperparameter selection. If params does not have a field called `num_inner_folds`, this script uses the same number of outer-loop folds as the number of inner-loop folds.

   - `feature_mat` (compulsory)
     A cell array of the full paths of the feature files storing a a cell array. For each feature file, a matrix "feature_mat" is assumed to be saved in this file. "feature_mat" can be a 2-D matrix with dimension of #features x #subjects or a 3-D matrix with dimension of #ROIs1 x #ROIs2 x #subjects. If "feature_mat" is 3-D, it is a connectivity matrix between two sets of ROIs, and only the lower-triangular off-diagonal entries will be considered as features because the connectivity matrix is symmetric.

   - `covariates` (compulsory)
     A matrix of the covariates to be regressed out from y (dimension: #subjects x #regressors). If the users only want to demean, an empty matrix should be passed in; if the users do not want to regress anything, set `params.covariates` = 'NONE'.

   - `y` (compulsory)
     A matrix of the target variable to be predicted (dimension: #subjects x #behaviors).

   - `outdir` (compulsory)
     A string of the full-path output directory.
 
   - `outstem` (compulsory)
     A string appended to the filename to specify the output files (y after regression, accuracy files, ...). For example, if `outstem = 'AngAffect_Unadj'`, then the final optimal accuracy file will be names as `fullfile(outdir, 'results', 'optimal_acc', 'AngAffect_Unadj.mat')`, and the optimal model hyperparameters of k-th fold will be saved in `fullfile(outdir, 'params', ['dis_' k '_cv'], 'selected_parameters_AngAffect_Unadj.mat')`.

   - `with_bias` (optional)
     A string (choose from '0' or '1') or a scalar (choose from 0 or 1). 1 indicates that a constant bias for every subject should be estimated. 1 by default.

   - `ker_param` (optional)
     A mat file containing Kx1 structs with two fields each: type (either 'corr', 'Gaussian', 'Exponential')  and scale. K denotes the number of kernels.

   - `threshold` (optional)
     A string or scalar specifying the threshold used to binarize the predicted score when the original y is binary. A scalar "threshold" is assumed to be saved in this file. The default is 0.5.

   - `group_kernel` (optional)
     A #kernel_groups cell array. Each element in the cell array contains a vector of size 1 x #kernels_in_the_group. For example if 4 kernels are to be used to train the multi kernel model, but they can be grouped such that 3 of the kernels uses the same hyperparameter lambda for regularization, then I can pass in the 'group_kernel' as {[1,3,4],[2]}. This means that kernels 1, 3 and 4 are regularized using the same lambda while kernel 2 is regularized using a different lambda. Do note that kernel 1 would correspond to the first feature mat in 'feature_file' so on and so forth. If this file is not passed in and also does NOT exist in "setup_file" or group_kernel is empty ,i.e. '', then it is assumed that each kernel will be regularised independently.  

   - `gpso_dir` (optional)
     A string, the full path of the directory to store the cloned git repository for Gaussian process code. Current script needs the Gaussian-Process Surrogate Optimisation ([GPSO](https://github.com/jhadida/gpso)) package to optimize the objective function. This git repository and its dependency "deck" (https://github.com/jhadida/deck) need to be both cloned into the same folder, which is passed into current script through params.gpso_dir. Default is `$CBIG_CODE_DIR/external_packages/matlab/non_default_packages/Gaussian_Process`

   - `domain` (optional)
     A 1 x 2 vector specifying the search domain for the regularization hyperparameters of the multi KRR model.The first element of the vector stores the lower bound and the second element stores the upper bound of the search domain. The default is set to [0 20].

   - `metric` (optional)
     A string, specifying how the value of lambda will be optimised with. Can choose the following: `corr`, `COD`, `predictive_COD`, `MAE`, `MAE_norm`,`MSE`,`MSE_norm`. Metric is `predictive_COD` by default.

### Saved results

When `CBIG_MultiKRR_workflow.m` is run, the optimal test prediction accuracies ('acc') would be stored in the folder `[outdir '/final_result_' outstem '.mat']`. acc would be a cell array of size K where K is the number of kernel types (see ker_param description). Each cell array would contain a #test_fold by #target_behaviors matrix of optimal test prediction accuracies.

Additionally, the field `optimal_stats` is saved in `[outdir '/final_result_' outstem '.mat']`. This contains a cell array. Each element in the cell array stores a data structure. There are 2 fields in the data structure, 'value' and 'description'. 'description' describes which accuracy metric the 'value' belongs to. 'value' will be a # folds x #behavioral_measures vector storing the values of the optimal test prediction accuracies corresponding to the specified acuracy metric. Within each field, there will be a total number of 7 entries, each entry corresponding to the possible acc_metrics.

----
## Updates

- Release v0.22.1 (08/03/2022): Fix bug in regressing covariates in functions calling `CBIG_regress_X_from_y_test.m`.
- Release v0.17.3 (13/01/2021): Initial release of multiKRR package

----
## Bugs and questions
Please contact Jianzhong Chen at chenjianzhong1996@gmail.com or Leon Ooi at leonooiqr@gmail.com.
