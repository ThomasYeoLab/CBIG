# Kernel Ridge Regression

## References
+ Li J, Kong R, Liegeois R, Orban C, Tan Y, Sun N, Holmes AJ, Sabuncu MR, Ge T, Yeo BTT, [Global signal regression strengthens association between resting-state functional connectivity and behavior](https://doi.org/10.1016/j.neuroimage.2019.04.016), Neuroimage, 2019, 196:126-141
+ Kong R, Li J, Orban C, et al. [Spatial Topography of Individual-Specific Cortical Networks Predicts Human Cognition, Personality, and Emotion](https://academic.oup.com/cercor/article/29/6/2533/5033556). Cerebral Cortex, 29(6):2533-2551, 2019

---

**Kernel ridge regression (KRR)** is a machine learning algorithm that can be used to predict behavioral values (`y`) of a test subject based on the similarity of some features (`feature`) between this test subject and all training subjects. For example, in the case of fluid intelligence prediction (i.e. `y` = fluid intelligence), if the test subject's functional connectivity (FC; `feature` = FC) is very similar to a training subject's FC, then we believe the test subject's fluid intelligence score is very similar to the training subject. 

To reduce overfitting, an L2-regularization term will be included. Meanwhile, our code is able to support multiple types of kernel (Gaussian kernel, exponential kernel, or correlation kernel) with different scaling factors. The regularization hyperparameter and the kernel hyperparameters will be determined by validation procedures, which can be performed in two different ways:

1. K-fold cross-validation: K-fold cross-validation (K-fold CV) stream
2. Splitting the data into training, validation, and test set: Training, validation, and test stream


## K-fold Cross-validation (K-fold CV) stream

For each behavioral phenotype (or other target measure you want to predict), the subjects are split into K folds. Care should be taken so that family members are not split between folds. For each test fold, `num_inner_fold`-fold cross-validation is repeatedly applied to the remaining K-1 folds with different regularization and kernel hyperparameters (i.e., inner-loop cross-validation). In practice, we included an L2-regularization term to reduce overfitting. The L2-regularization parameter (as well as other hyperparameters) is determined via the inner-loop cross-validation procedure. The optimal hyperparameters from the inner-loop crossvalidation are then used to predict the behavioral phenotype in the test fold. Accuracy is measured by correlating the predicted and actual behavioral measure across all subjects within the test fold. By repeating the procedure for each test fold, each behavior yielded K correlation accuracies, which were then averaged across the K folds. Because a single K-fold cross-validation might be sensitive to the particular split of the data into folds, the above K-fold cross-validation can be **repeated multiple times** by passing in different fold-split file (see `sub_fold` below). For example, if we do 20-fold cross-validation, the process can be repeated with 100 different splits.

Steps of this stream:
* Regress nuisance variables (e.g. head movement) from the target variable (`y`)
* Generate the kernels for inner-loop CV step and the test CV step
* Inner-loop CV with all possible hyperparameter combinations
* Test CV with the same sets of hyperparameter
* Based on the test accuracies in inner-loop CV, select the optimal hyperparameter setting and obtain the corresponding test accuracy in the test CV. 

### How to use the scripts

There are two sets of scripts you can use. The top-level wrapper functions of each set are `CBIG_KRR_workflow.m` and `CBIG_KRR_workflow_LITE.m` respectively. The usages of these two wrapper functions are the same. Therefore, the following sections will only discuss about `CBIG_KRR_workflow.m`. The functionality of these two sets of scripts (with or without `LITE` appendix) are slightly different. When `CBIG_KRR_workflow.m` is used, the functional similarity matrices for each fold are saved. If the user chooses Gaussian or exponential kernel, the mean and standard deviation of features are computed from the training subjects and applied on the test subjects for each fold. However, when `CBIG_KRR_workflow_LITE.m` is used, the functional similarity matrix is only computed and saved once across all subjects so that less disk space is needed. In the following training and testing on each fold, the scripts grab the similarity matrix of current fold by indexing the entire similarity matrix. **`CBIG_KRR_workflow_LITE.m` is only able to support correlation kernel**.

To use `CBIG_KRR_workflow.m`, you can either pass in a single setup structure or file (the first argument `setup_param` of this function), or pass in a set of individual parameters (Compulsory variables and optional variables via varargin).

For example, 
1) If we wish to just pass in setup_param 
CBIG_KRR_workflow( setup_param)


2) If we wish to pass in input variables one by one and also specifying optional variables `with_bias` and `lambda_set_file`
CBIG_KRR_workflow( '', 1, sub_fold_file, y_file, covariate_file, feature_file, num_inner_folds, outdir, outstem, 'with_bias', 0, 'lambda_set_file', 'xxx/xxx/lambda_set_file.mat')

#### 1. `setup_param`
Here is an example file for setup_param: `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples/output/KernelRidgeRegression/setup_file.mat`

If the user decided to pass in the setup structure or file, then input variables from `sub_fold file` to `varargin` are not needed. `setup_param` can be setup parameter (a Matlab structure variable) or full path of the setup parameter (.mat storing a Matlab structure variable called `params`). The structure should have the following structure fields:
* `param.sub_fold`
  * `param.sub_fold` is a structure which tells how the subjects are split. 
  * `param.sub_fold(i).fold_index` is a N x 1 binary vector, where N is the total number of subjects. `param.sub_fold(i).fold_index(j) = 1` indicates that the `j`-th subject is in the test set of the `i`-th outer-loop CV iteration, while `param.sub_fold(i).fold_index(j) = 0` indicates that the `j`-th subject is in the training set of the `i`-th outer-loop CV iteration.
  * The user can use `$CBIG_CODE_DIR/utilities/matlab/predictive_models/utilities/CBIG_cross_validation_data_split.m` to generate the structure variable `sub_fold`. Here is an example of using this function. 
    
    Inputs:
    * `subject_list`: a text file containing all subject IDs (each line corresponds to one subject ID). 
    * `family_csv`: a CSV file where there is a column representing the family ID of each subject. If all the subjects are unrelated, you can pass in `'none'`. 
    * `subject_header`: the header of the subject ID column in `family_csv`. If `family_csv = 'none'`, you can pass in an empty string `''`.
    * `family_header`: the header of the family ID column in `family_csv`. If `family_csv = 'none'`, you can pass in an empty string `''`. 
    * `num_folds`: how many folds that the subjects will be split into. 
    * `seed`: a random seed used to initialize the random seed generator. You can use different random seeds to split the subjects differently. 
    * `outdir`: a user-defined output directory. A file containing the `sub_fold` variable will be saved as `[outdir '/no_relative_' num_folds '_fold_sub_list.mat]`.
    * `delimiter`: the delimiter used in `family_csv`.
    
    Commands:
    ```
    subject_list = '<your_own_path>/subject_list.txt';
    family_csv = '<your_own_path>/familyID.csv';
    subject_header = 'Subject_ID';
    family_header = 'family_ID';
    num_folds = 20;
    seed = 1;
    outdir = '<your_own_path>/outdir';
    delimiter = ',';
    sub_fold = CBIG_cross_validation_data_split( subject_list, family_csv, subject_header, family_header, num_folds, seed, outdir, delimiter );
    ``` 
* `param.y`
  * `param.y` is a matrix of behavioral measures (or other measures you want to predict) with size N x M, where N is the total number of subjects, and M is the number of behavioral measures. 
  * If your behavioral measures are stored in several CSV files, we provide a function `$CBIG_CODE_DIR/utilities/matlab/predictive_models/utilities/CBIG_read_y_from_csv.m` to extract the behavioral measures from these CSV files. Here is an example of using this function.
    
    Inputs:
    * `csv_files`: the filenames of the CSV files containing the behavioral measures.
    * `subject_header`: the header of the subject ID column in the CSV files. It can change according to the header information in the CSV files.
    * `y_names`: the name of the behavioral measures you want to predict.
    * `y_types`: determines how the behavioral scores will be read in. `y_types{i} = 'continuous'` assumes the `i`-th behavioral scores listed in the CSV files are continuous values. Sometimes, there could be categorical measures written as strings in the CSV file (e.g. Alzhemier's disease diagnosis could be labeled as 'AD', 'MCI' or 'control'). In this case, you can set `y_types{i} = 'categorical'` for this measure. The number of cell arrays in `y_types` should match the number of cell arrays in `y_names`.
    * `subject_list`: a text file containing all subject IDs (each line corresponds to one subject ID). 
    * `outname` is the output filename. The matrix `y` will be saved in this file.
    * `delimiter`: pass in a string as the delimiter if the columns in all the CSV files can be separated by the same `delimiter`. If the delimiters in the CSV files are different, you can use cell arrays to pass in the delimiters of each CSV file, e.g. `{',', ';'}` means the first CSV file is delimited by `,` and the second CSV file is delimited by `;`, and so on and so forth.
    
    Commands:
    ```
    csv_files = {'<your_own_path>/CSV1.csv', '<your_own_path>/CSV2.csv'};
    subject_header = 'Subject_ID';
    y_names = {'Behavior_1', 'Behavior_2'};
    y_types = {'continuous', 'categorical'};
    subject_list = '<your_own_path>/subject_list.txt';
    outname = '<your_own_path>/outname.mat';
    delimiter = ','
    y = CBIG_read_y_from_csv( csv_files, subject_header, y_names, y_types, subject_list, outname, delimiter );
    ```
* `param.covariates`
  * `param.covariates` is a matrix of the covariates you want to regress out from the behavioral measures. The size is N x P, where N is the total number of subjects, and P is the number of covariate variables you want to regress. (Note that the number of regressors can be different from the number of covariate measures, if there is categorical covariate measure. See the description of `covariate_types` below.)
  * If your covariates are stored in several CSV files, we provide a function `$CBIG_CODE_DIR/utilities/matlab/predictive_models/utilities/CBIG_generate_covariates_from_csv.m` to extract the covariate measures from these CSV files. Here is an example of using this function.
    
    Inputs:
    * `csv_files`: the filenames of the CSV files containing the covariate measures.
    * `subject_header`: the header of the subject ID column in the CSV files. It can change according to the header information in the CSV files.
    * `covariate_names`: the name of the covariate measures you want to regress.
    * `covariate_types`: the type of  covariate values that will be read in. `covariate_types{p} = 'continuous'` assumes the `p`-th covariate values are continuous numbers in the CSV files. You might also want to regress categorical covariates (e.g. race could be labeled as 'Asian', 'Hispanic', or 'African'). In this case, you can set `covariate_types{p} = 'categorical'` for this measure. An N x Q_p regressor matrix will be generated for this covariate measure, where Q_p is the total number of categories of the `p`-th covariate minus 1. The number of cell arrays in `covariate_types` should equal to the number of cell arrays in `covariate_names`.
    * `subject_list`: a text file containing all subject IDs (each line corresponds to one subject ID). 
    * `FD_file`: If there is a need to regress the framewise displacement (FD), i.e. `covariate_names` contains a cell array called `'FD'`, `FD_file` is the name of a text file in which each line corresponds to the mean FD value of a subject. If `covariate_names` does not contain `'FD'`, `FD_file` will be ignored so that you can pass in `''` or `'none'`.
    * `DVARS_file`: If there is a need to regress the DVARS of each subject, i.e. `covariate_names` contains a cell array called `'DVARS'`, `DVARS_file` is the name of a text file in which each line corresponds to the mean DVARS value of a subject. If `covariate_names` does not contain `'DVARS'`£¬ `DVARS_file` will be ignored so that you can pass in `''` or `'none'`.
    * `outname` is the output filename. The matrix `covariates` will be saved in this file.
    * `delimiter`: pass in a string as the delimeter if the columns in all the CSV files can be separated by the same `delimiter`. If the delimiters in the CSV files are different, you can use cell arrays to pass in the delimiters of each CSV file, e.g. `{',', ';'}` means the first CSV file is delimited by `,` and the second CSV file is delimited by `;`, and so on and so forth.
    
    Commands:
    ```
    csv_files = {'<your_own_path>/CSV1.csv', '<your_own_path>/CSV2.csv'};
    subject_header = 'Subject_ID';
    covariate_names = {'Age', 'Sex', 'Race'};
    covariate_types = {'continuous', 'categorical', 'categorical'};
    subject_list = '<your_own_path>/subject_list.txt';
    FD_file = '<your_own_path>/FD_file.txt';
    DVARS_file = '<your_own_path>/DVARS_file.txt';
    outname = '<your_own_path>/outname.mat';
    delimiter = ',';
    covariates = CBIG_generate_covariates_from_csv( csv_files, subject_header, covariate_names, covariate_types, subject_list, FD_file, DVARS_file, outname, delimiter );
    ```
* `param.feature_mat`
  * It is a matrix of the features used as the explanatory variables in the prediction.
  * `feature_mat` can be a 2-D matrix with dimension of F x N, where F is the number of features and N is the number of subjects. In this case, each column is a feature vector of a single subject.
  * `feature_mat` can also be a 3-D matrix with dimension of R1 x R2 x N, where R1 is the number of ROIs_1, R2 is the number of ROIs_2, and N is the number of subjects. In this case, it is a connectivity matrix between ROIs_1 and ROIs_2. Only the lower-triangular off-diagonal entries will be considered.
* `param.num_inner_folds`
  * A scalar, the number of inner-loop cross-validation folds. For example, `20`.
* `param.outdir`
  * It is a string of the output directory (full path). An example folder structure in the output directory can be found in `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples/output/KernelRidgeRegression/`.
* `param.outstem`
  * It is a string appended to all the output filenames. For example, you may want to run this workflow multiple times for the same dataset but with different target variables, different covariates, ... `outstem` can help you to differentiate these multiple runs even if they are saved in the same output folder. For instance, the filename of the final test accuracies with optimal hyperparameters will be `[param.outdir '/final_result_' param.outstem '.mat']` if `param.outstem` is not empty. If `param.outstem` is empty, the filename will become `[param.outdir '/final_result.mat']`.
* `param.with_bias`
  * It is a string or a scalar that can be either 0 or 1. Default is 1. If `param.with_bias = 1`, the cost function to be minimized is `(y - K*alpha - beta)^2 + (regularization of alpha)`, where `beta` is a constant bias for every subject, estimated from the data. If `param.with_bias = 0`, the cost function to be minimized is `(y - K*alpha)^2 + (regularization of alpha)`. It means that the kernel regression fitted curve must pass through the origin since the bias term `beta` (acting as an intercept) is not included. If your target measures are not demeaned before and you still set `param.with_bias = 0`, your predicted scores would not be accurate. Hence `param.with_bias = 0` is not recommended.
* `param.save_kernel`
  * **Only for `LITE` version** It is a string or a scalar that can be either 0 or 1. Default is 1. `save_kernel = 0` (or '0') means the algorithm is do not save kernel into files. `save_kernel = 1` (or '1') means the algorithm is save kernel into files.
* `param.ker_param`
  * It is a structure of length L, where L is the number of kernel parameters the user would like to pass into the workflow.
  * It should contain two subfields: `type` and `scale`.
  * `ker_param(l).type` is a string of the type of the l-th kernel. Options include: 
    * `'corr'`             - Pearson's correlation;
    * `'Gaussian'`         - Gaussian kernel;
    * `'Exponential'`      - exponential kernel.
  * `ker_param(l).scale` is a scalar specifying the scaling factor of the l-th kernel (only applicable for Gaussian kernel and exponential kernel). If `ker_param(l).type == 'corr'`, then set `ker_param(l).scale = NaN`.
* `param.lambda_set`
  * It is a vector of numbers for grid search of `lambda` (the regularization parameter). For example, `[0.001 0.005 0.01:0.02:0.1 0.2:0.1:1 5:5:20]`.
* `param.threshold_set`
  * The target measure you want to predict can be binary (e.g. sex). In this case, we introduce another hyperparameter, `threshold`, as a separation point to predict binary measures.
  * `threshold_set` is a vector of numbers (between -1 and 1) used for grid search of `threshold`.
* `param.metric`
  * A string indicating the metric used to define prediction loss. The loss is used to choose hyperparameters.Options include: 
    * `'corr'`             - Pearson's correlation;
    * `'COD'`              - Coefficient of determination. Defined as 1-||y_pred-y_test||^2/||mean(y_test)-y_test||^2, where y_pred is the prediction of the test data, y_test is the groud truth of the test data, and mean(y_test) is the mean of test data.
    * `'Predictive_COD'`   - Predictive coefficient of determination. Defined as 1-||y_pred-y_test||^2/||mean(y_train)-y_test||^2, where y_pred is the prediction of the test data, y_test is the groud truth of the test data, and mean(y_train) is the mean of training data.
    * `'MAE'`              - mean absolute error.
    * `'MAE_norm'`         - mean absolute error divided by the standard derivation of the target variable of the training set.
    * `'MSE'`              - mean squared error.
    * `'MSE_norm'`         - mean squared error divided by the variance of the target variable of the training set.



#### 2. Compulsory and Optional Variables if `setup_param` is not passed in
If the user did not prepare the setup file, he/she needs to pass in the parameters needed one by one through the compulsory and optional variables (after the input variable `save_setup`) and leave `setup_param` as empty. The first 7 Variables are compulsory (ie from `sub_fold_file` to `outstem`). Subsequently `varargin` will parse the optional input variables. In Matlab, `varargin` grabs all the parameters passed in, starting from the position of `varargin` till the last input, and store them in cell arrays. In the case of `CBIG_KRR_workflow.m`, since `varargin` is the 10th argument, the 10th to the last inputs will all be stored in `varargin`. If the user wishes to specify certain optional variables like 'metric', all the user needs to do is (...,outstem, 'metric', metric). The user needs to specify the variable name first ('metric') and then followed by the variable value (metric).

Overall, the parameters you need to pass in are:
* sub_fold_file
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a structure variable called `sub_fold`.
  * `sub_fold` is a structure variable which tells how the subjects are split. 
  * `sub_fold(i).fold_index` is a N x 1 binary vector, where N is the total number of subjects. `sub_fold(i).fold_index(j) = 1` indicates that the `j`-th subject is in the test set of the `i`-th outer-loop CV iteration, while `sub_fold(i).fold_index(j) = 0` indicates that the `j`-th subject is in the training set of the `i`-th outer-loop CV iteration.
  * The user can use `$CBIG_CODE_DIR/utilities/matlab/predictive_models/utilities/CBIG_cross_validation_data_split.m` to generate the structure variable `sub_fold`. Here is an example of using this function. 
    
    Inputs:
    * `subject_list`: a text file containing all subject IDs (each line corresponds to one subject ID). 
    * `family_csv`: a CSV file in which there is a column representing the family ID of each subject. If all the subjects are unrelated, you can pass in `'none'`. 
    * `subject_header`: the header of the subject ID column in `family_csv`. If `family_csv = 'none'`, you can pass in an empty string `''`.
    * `family_header`: the header of the family ID column in `family_csv`. If `family_csv = 'none'`, you can pass in an empty string `''`. 
    * `num_folds`: how many folds that the subjects will be split into. 
    * `seed`: a random seed used to initialize the random seed generator. You can use different random seeds to split the subjects differently. 
    * `outdir`: a user-defined output directory. A file containing the `sub_fold` variable will be saved as `[outdir '/no_relative_' num_folds '_fold_sub_list.mat]`.
    * `delimiter`: the delimiter used in `family_csv`.
    
    Commands:
    ```
    subject_list = '<your_own_path>/subject_list.txt';
    family_csv = '<your_own_path>/familyID.csv';
    subject_header = 'Subject_ID';
    family_header = 'family_ID';
    num_folds = 20;
    seed = 1;
    outdir = '<your_own_path>/outdir';
    delimiter = ',';
    sub_fold = CBIG_cross_validation_data_split( subject_list, family_csv, subject_header, family_header, num_folds, seed, outdir, delimiter );
    ``` 
* y_file
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a matrix called `y`.
  * `y` is a matrix of behavioral measures (or other measures you want to predict) with size N x M, where N is the total number of subjects, and M is the number of behavioral measures. 
  * If your behavioral measures are stored in several CSV files, we provide a function `$CBIG_CODE_DIR/utilities/matlab/predictive_models/utilities/CBIG_read_y_from_csv.m` to extract the behavioral measures from these CSV files. Here is an example of using this function.
    
    Inputs:
    * `csv_files`: the filenames of the CSV files containing the behavioral measures.
    * `subject_header`: the header of the subject ID column in the CSV files. It can change according to the header information in the CSV files.
    * `y_names`: the name of the behavioral measures you want to predict.
    * `y_types`: determines how the behavioral scores will be read in. `y_types{i} = 'continuous'` assumes the `i`-th behavioral scores listed in the CSV files are continous values. Sometimes, there could be categorical measures written as strings in the CSV file (e.g. Alzhemier's disease diagnosis could be labeled as 'AD', 'MCI' or 'control'). In this case, you can set `y_types{i} = 'categorical'` for this measure. The number of cell arrays in `y_types` should match the number of cell arrays in `y_names`.
    * `subject_list`: a text file containing all subject IDs (each line corresponds to one subject ID). 
    * `outname` is the output filename. The matrix `y` will be saved in this file.
    * `delimiter`: pass in a string as the delimiter if the columns in all the CSV files are separated by the same `delimiter`. If the delimiters in the CSV files are different, you can use cell arrays to pass in the delimiters of each CSV file, e.g. `{',', ';'}` means the first CSV file is delimited by `,` and the second CSV file is delimited by `;`, and so on and so forth.
    
    Commands:
    ```
    csv_files = {'<your_own_path>/CSV1.csv', '<your_own_path>/CSV2.csv'};
    subject_header = 'Subject_ID';
    y_names = {'Behavior_1', 'Behavior_2'};
    y_types = {'continuous', 'categorical'};
    subject_list = '<your_own_path>/subject_list.txt';
    outname = '<your_own_path>/outname.mat';
    delimiter = ','
    y = CBIG_read_y_from_csv( csv_files, subject_header, y_names, y_types, subject_list, outname, delimiter );
    ```
* covariate_file
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a matrix called `covariates`.
  * `covariates` is a matrix of the covariates you want to regress out from the behavioral measures. The size is N x P, where N is the total number of subjects, and P is the number of covariate variables you want to regress. (Note that the number of regressors can be different from the number of covariate measures, if there is categorical covariate measure. See the description of `covariate_types` below.)
  * If your covariates are stored in several CSV files, we provide a function `$CBIG_CODE_DIR/utilities/matlab/predictive_models/utilities/CBIG_generate_covariates_from_csv.m` to extract the covariate measures from these CSV files. Here is an example of using this function.
    
    Inputs:
    * `csv_files`: the filenames of the CSV files containing the covariate measures.
    * `subject_header`: the header of the subject ID column in the CSV files. It can change according to the header information in the CSV files.
    * `covariate_names`: the name of the covariate measures you want to regress.
    * `covariate_types`: the type of covariates values that will be read in. `covariate_types{p} = 'continuous'` assumes the `p`-th covariate values are continuous numbers in the CSV files. You might also want to regress categorical covariates (e.g. race could be labeled as 'Asian', 'Hispanic', or 'African'). In this case, you can set `covariate_types{p} = 'categorical'` for this measure. An N x Q_p regressor matrix will be generated for this covariate measure, where Q_p is the total number of categories of the `p`-th covariate minus 1. The number of cell arrays in `covariate_types` should equal to the number of cell arrays in `covariate_names`.
    * `subject_list`: a text file containing all subject IDs (each line corresponds to one subject ID). 
    * `FD_file`: If there is a need to regress the framewise displacement (FD), i.e. `covariate_names` contains a cell array called `'FD'`, `FD_file` is the name of a text file in which each line corresponds to the mean FD value of a subject. If `covariate_names` does not contain `'FD'`, `FD_file` will be ignored so that you can pass in `''` or `'none'`.
    * `DVARS_file`: If there is a need to regress the DVARS of each subject, i.e. `covariate_names` contains a cell array called `'DVARS'`, `DVARS_file` is the name of a text file in which each line corresponds to the mean DVARS value of a subject. If `covariate_names` does not contain `'DVARS'`. `DVARS_file` will be ignored so that you can pass in `''` or `'none'`.
    * `outname` is the output filename. The matrix `covariates` will be saved in this file.
    * `delimiter`: pass in a string as the delimeter if the columns in all the CSV files can be separated by the same `delimiter`. If the delimiters in the CSV files are different, you can use cell arrays to pass in the delimiters of each CSV file, e.g. `{',', ';'}` means the first CSV file is delimited by `,` and the second CSV file is delimited by `;`, and so on and so forth.
    
    Commands:
    ```
    csv_files = {'<your_own_path>/CSV1.csv', '<your_own_path>/CSV2.csv'};
    subject_header = 'Subject_ID';
    covariate_names = {'Age', 'Sex', 'Race'};
    covariate_types = {'continuous', 'categorical', 'categorical'};
    subject_list = '<your_own_path>/subject_list.txt';
    FD_file = '<your_own_path>/FD_file.txt';
    DVARS_file = '<your_own_path>/DVARS_file.txt';
    outname = '<your_own_path>/outname.mat';
    delimiter = ',';
    covariates = CBIG_generate_covariates_from_csv( csv_files, subject_header, covariate_names, covariate_types, subject_list, FD_file, DVARS_file, outname, delimiter );
    ```
* feature_file
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a matrix called `feature_mat`.
  * `feature_mat` is a matrix of the features used as the explanatory variables in the prediction.
  * `feature_mat` can be a 2-D matrix with dimension of F x N, where F is the number of features and N is the number of subjects. In this case, each column is a feature vector of a single subject.
  * `feature_mat` can also be a 3-D matrix with dimension of R1 x R2 x N, where R1 is the number of ROIs_1, R2 is the number of ROIs_2, and N is the number of subjects. In this case, it is a connectivity matrix between ROIs_1 and ROIs_2. Only the lower-triangular off-diagonal entries will be considered.
* num_inner_folds
  * A scalar, the number of inner-loop cross-validation folds. For example, `20`.
* outdir
  * It is a string of the output directory (full path). An example folder structure in the output directory can be found in `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples/output/KernelRidgeRegression/`.
* outstem
  * It is a string appended to all the output filenames. For example, you may want to run this workflow multiple times for the same dataset but with different target variables, different covariates, ... This can help you to differentiate these multiple runs even if they are saved in the same output folder. For instance, the filename of the final test accuracies with optimal hyperparameters will be `[outdir '/final_result_' outstem '.mat']` if `outstem` is not empty. If `outstem` is empty, the filename will become `[outdir '/final_result.mat']`.
* (...,'with_bias', with_bias,...)
  * It is a string or a scalar that can be either 0 or 1. Default is 1. If `with_bias = 1`, the cost function to be minimized is `(y - K*alpha - beta)^2 + (regularization of alpha)`, where `beta` is a constant bias for every subject, estimated from the data. If `with_bias = 0`, the cost function to be minimized is `(y - K*alpha)^2 + (regularization of alpha)`. It means that the kernel regression fitted curve must pass through the origin since the bias term `beta` (acting as an intercept) is not included. If your target measures are not demeaned before and you still set `with_bias = 0`, your predicted scores would not be accurate. Hence `with_bias = 0` is not recommended.
* (...,'save_kernel', save_kernel,...)
  * **Only for `LITE` version** It is a string or a scalar that can be either 0 or 1. Default is 1. `save_kernel = 0` (or '0') means the algorithm is do not save kernel into files. `save_kernel = 1` (or '1') means the algorithm is save kernel into files.
* (...,'ker_param_file', ker_param_file,...)
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a structure called `ker_param`.
  * `ker_param` is a structure of length L, where L is the number of kernel parameters the user would like to pass into the workflow.
  * `ker_param` should contain two subfields: `type` and `scale`.
  * `ker_param(l).type` is a string of the type of the l-th kernel. Options include: 
    * `'corr'`             - Pearson's correlation;
    * `'Gaussian'`         - Gaussian kernel;
    * `'Exponential'`      - exponential kernel.
  * `ker_param(l).scale` is a scalar specifying the scaling factor of the l-th kernel (only applicable for Gaussian kernel and exponential kernel). If `ker_param(l).type == 'corr'`, then set `ker_param(l).scale = NaN`.
  * If ker_param_file is not passed in, the default `ker_param` will be set as `ker_param.type = 'corr'; ker_param.scale = NaN;`.
* (...,'lambda_set_file', lambda_set_file,...)
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a vector called `lambda_set`.
  * `lambda_set` is a vector of numbers for grid search of `lambda` (the regularization parameter). For example, `[0.001 0.005 0.01:0.02:0.1 0.2:0.1:1 5:5:20]`.
  * If `lambda_set_file` is not passed in, the default `lambda_set` will be set as `lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20]`.
* (...,'threshold_set_file', threshold_set_file,...)
  * A string, which is the full-path name of a `.mat` file. The `.mat` file contains a vector called `threshold_set`.
  * The target measure you want to predict can be binary (e.g. sex). In this case, we introduce another hyperparameter, `threshold`, as a separation point to predict binary measures.
  * `threshold_set` is a vector of numbers (between -1 and 1) used for grid search of `threshold`.
  * If `threshold_set_file` is not passed in, the default `threshold_set` will be set as `threshold_set = [-1:0.1:1]`.
* (...,'metric',metric, ...)
  * A string indicating the metric used to define prediction loss. The loss is used to choose hyperparameters.Options include: 
    * `'corr'`             - Pearson's correlation;
    * `'COD'`              - Coefficient of determination. Defined as 1-||y_pred-y_test||^2/||mean(y_test)-y_test||^2, where y_pred is the prediction of the test data, y_test is the groud truth of the test data, and mean(y_test) is the mean of test data.
    * `'Predictive_COD'`   - Predictive coefficient of determination. Defined as 1-||y_pred-y_test||^2/||mean(y_train)-y_test||^2, where y_pred is the prediction of the test data, y_test is the groud truth of the test data, and mean(y_train) is the mean of training data.
    * `'MAE'`              - mean absolute error.
    * `'MAE_norm'`         - mean absolute error divided by the standard derivation of the target variable of the training set.
    * `'MSE'`              - mean squared error.
    * `'MSE_norm'`         - mean squared error divided by the variance of the target variable of the training set.
  * If `metric` is not passed in, the default `metric` will be set as `metric = 'predictive_COD'`.
## Training, validation, and test stream

If you have enough number of observations (e.g. hundreds of thousands of subjects), you might want to simply divide your dataset into training, validation and test sets, instead of K-fold cross-validation. In this case, the hyperparameters will be determined based on the optimal accuracies in the validation set. After that, the model parameters (i.e. the kernel regression coefficients) will be obtained from the training set with the optimal hyperparameters, and then applied to the test set to get the test accuracy.

Brief steps include:
* Regress the nuisance variables from the target variables `y`
* Generate kernels for training, validation, and test phases
* For each possible hyperparameter combination, train the model on the training set, and obtain the accuracies and predicted values in the validation set
* For each possible hyperparameter combination, train the model on the training set, and obtain the accuracies and predicted values in the test set
* Based on the accuracies in the validation set, select the optimal hyperparameters. Extract the accuracies and predicted values in the test set corresponding to the optimal hyperparameters.

### How to use the scripts

The usage of `CBIG_KRR_workflow.m` is roughly the same as the cross-validation stream. The only difference is how to specify `sub_fold`. In this case, the length of `sub_fold` is 1. `sub_fold.fold_index` is still a #subjects x 1 vector. However, the values of `sub_fold.fold_index` can be 0, 1 or 2:
```
0 - training set
1 - test set
2 - validation set
```
For the specification of other parameters, please refer to the cross-validation section.

### Saved results

When `CBIG_KRR_workflow.m` is run, a folder named `results` will be generated in the path specified in `outdir`. The `results` folder should contain 5 sub-folders: 

   - `y`
     This folder contains the original y values and the y values after the covariates have been regressed for each fold. 

   - `FSM_innerloop`
     This folder contains the kernels for training subjects only. It should be a #training subjects x #training subjects matrix. The kernels are saved in a separate folder for each fold.

   - `FSM_test`
     This folder contains the kernels for all subjects. It should be a #subjects x #subjects matrix. The kernels are saved in a separate folder for each fold.

   - `innerloop_cv`
     This folder contains the accuracy, loss and predicted y value for each lambda value for the innerloop. The results are saved in a separate mat file for each fold.

   - `test_cv`
     This folder contains the accuracy, loss and predicted y value for each lambda value for the outerloop. The results are saved in a separate mat file for each fold.

The final accuracies for the test folds are saved in `outdir/results` under `final_result_<outstem>.mat` with the following fields: 

   - `optimal_acc`
     Test accuracy for each fold (given in correlation).

   - `optimal_kernel`
     A #outerfolds x #behaviors struct containing the kernel type selected for each test fold.

   - `optimal_lambda`
     A #outerfolds x #behaviors matrix containing lambda selected for each test fold.

   - `optimal_threshold`
     TA #outerfolds x #behaviors matrix containing the threshold selected for each test fold.

   - `y_predict_concat` 
     A #subjects x #behaviors matrix of predicted target values.

   - `optimal_stats`
     A cell array storing the accuracies of each possible accuracy metric (eg. corr, MAE, etc). Each cell array is #outerfolds x #behaviors.

## Example usage

We provided a fake example to show how to use this workflow (for the K-fold cross-validation stream). See the readme file: `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples/README.md`.

The script can be found as `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples/scripts/CBIG_LiGSR_example_KRR.sh`

You can check the output folder structure in this directory: `$CBIG_CODE_DIR/stable_projects/preprocessing/Li2019_GSR/examples/output/KernelRidgeRegression`.

## Updates

- Release v0.22.6 (11/05/2022): Removed option of Gaussian and Exponential kernels from `CBIG_KRR_generate_kernels_LITE.m`.
- Release v0.22.1 (08/03/2022): Fix bug in regressing covariates in functions calling `CBIG_regress_X_from_y_test.m`.
- Release v0.21.2 (09/12/2021): Add functionality to regress covariates from features within cross-validation stream.
- Release v0.15.4 (05/11/2019): Change default lambda set; add hyperparameter-tuning metric.
- Release v0.15.3 (16/10/2019): Add references; add LITE version of kernel regression code.
- Release v0.14.1 (03/09/2019): Optimized the speed of kernel regression scripts.
- Release v0.13.1 (19/07/2019):
  1. Added bias term in kernel regression cost function.
  2. Added unit test script for general kernel regression package.
- Release v0.9.4 (01/04/2019): Renamed matlab variable `feature` to `feature_mat` to avoid clashing with matlab built-in function.
- Release v0.9.0 (13/02/2019): Initial release of general kernel regression package

## Bugs and questions
Please contact Jingwei Li at jingweili.sjtu.nus@gmail.com and Ru(by) Kong at roo.cone@gmail.com.