function CBIG_KRR_example_wrapper_input_args(output_dir)

% CBIG_KRR_example_input_args(output_dir)
% This wrapper function performs KRR using a publicly-available air quality dataset.
% The aim of this example to let the users familiarize our KRR workflow, please
% refrain from making any conclusion regarding this air quality dataset per se.
% This script will generate all necessary input arguments and pass them to KRR workflow.
% In this example, KRR is run using 5-fold cross validation with 3 different splits
% Input:
%  - output_dir:
%    the absolute path for output directory 
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% the absolute path for all input files
input_path = [getenv("CBIG_CODE_DIR") '/utilities/matlab/predictive_models/KernelRidgeRegression/example/input'];

% the script runs KRR with 3 different splits
for split = 1:3
    
    %outdir
    outdir_each_split = [output_dir '/' num2str(split)];
    
    %sub_fold
    subject_list = [input_path '/data_id.txt'];
    family_csv = 'none';
    subject_header = '';
    family_header = '';
    num_folds = 5;
    CBIG_cross_validation_data_split(subject_list, family_csv, subject_header, family_header, ...
        num_folds, split, outdir_each_split);
    sub_fold = [outdir_each_split '/no_relative_5_fold_sub_list.mat'];
    
    %y
    csv_files = {[input_path '/CO_groundtruth.csv']};
    subject_header = 'data_id';
    y_names = { 'CO(GT)'};
    y_types = {'continuous'};
    outname = [outdir_each_split '/CO_groundtruth.mat'];
    delimiter = ',';
    CBIG_read_y_from_csv( csv_files, subject_header, y_names, y_types, subject_list, outname, delimiter);
    y = [outdir_each_split '/CO_groundtruth.mat'];

    %covariates
    csv_files = {[input_path '/date.csv']};
    subject_header = 'data_id';
    covariate_names = {'Date'};
    covariate_types = {'continuous'};
    FD_file = 'none';
    DVARS_file = 'none';
    outname = [outdir_each_split '/air_quality_covariate.mat'];
    delimiter = ',';
    CBIG_generate_covariates_from_csv( csv_files, subject_header, covariate_names, ...
       covariate_types, subject_list, FD_file, DVARS_file, outname, delimiter );
    covariates = [outdir_each_split '/air_quality_covariate.mat'];

    %feature_mat
    feature_mat = [input_path '/feature.mat'];  
    
    %inner folds
    num_inner_folds = 5;
    
    %outstem
    outstem = ['CO_' num2str(split)];   

    %metric
    metric = 'predictive_COD';
    
    CBIG_KRR_workflow_LITE('', 1, sub_fold, y, covariates, feature_mat, num_inner_folds, ...
        outdir_each_split, outstem, 'with_bias', 1, 'metric', metric);
end

% check result
CBIG_KRR_example_check_result(output_dir);   

end
