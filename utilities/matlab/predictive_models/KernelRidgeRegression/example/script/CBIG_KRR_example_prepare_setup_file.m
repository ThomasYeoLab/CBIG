function CBIG_KRR_example_prepare_setup_file(input_path, outdir, split)

% CBIG_KRR_example_prepare_setup_file(input_path, outdir, split)
% This function prepares a setup file for KRR workflow. The generated setup file is
% only applicable for the purpose of this KRR example.
% Input:
%  - input_path:
%    the absolute path of the input directory
%  - outdir:
%    the absolute path for output directory
%  - split:
%    a number that is used to initialise the random number generator to split the dataset
%    note that the variable type should be 'string'
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%param.outdir
outdir_each_split = [outdir '/' split];
param.outdir = outdir_each_split;

%param.sub_fold
subject_list = [input_path '/data_id.txt'];
family_csv = 'none';
subject_header = '';
family_header = '';
num_folds = 5;
sub_fold = CBIG_cross_validation_data_split(subject_list, family_csv, subject_header, family_header, ...
    num_folds, str2num(split), outdir_each_split);
param.sub_fold = sub_fold;

%param.y
csv_files = {[input_path '/CO_groundtruth.csv']};
subject_header = 'data_id';
y_names = { 'CO(GT)'};
y_types = {'continuous'};
outname = [outdir_each_split '/CO_groundtruth.mat'];
delimiter = ',';
y = CBIG_read_y_from_csv( csv_files, subject_header, y_names, y_types, subject_list, outname, delimiter);
param.y = y;

%param.covariates
csv_files = {[input_path '/date.csv']};
subject_header = 'data_id';
covariate_names = {'Date'};
covariate_types = {'continuous'};
FD_file = 'none';
DVARS_file = 'none';
outname = [outdir_each_split '/air_quality_covariate.mat'];
delimiter = ',';
covariates = CBIG_generate_covariates_from_csv( csv_files, subject_header, covariate_names, ...
   covariate_types, subject_list, FD_file, DVARS_file, outname, delimiter );
param.covariates = covariates;

%param.feature_mat
feature = load([input_path '/feature.mat']);
param.feature_mat = feature.feature_mat;

%param.num_inner_folds
param.num_inner_folds = 5;

%param.outstem
param.outstem = ['CO_' split];

%param.with_bias
param.with_bias = 1;

%param.ker_param
param.ker_param.type = 'corr';
param.ker_param.scale = NaN;

%param.lambda_set
param.lambda_set = [0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];

%param.metric
param.metric = 'predictive_COD';

param_name = [outdir_each_split '/setup_' split '.mat'];
save (param_name, '-struct', 'param');
