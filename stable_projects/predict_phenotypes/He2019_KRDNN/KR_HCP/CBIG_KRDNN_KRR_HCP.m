function CBIG_KRDNN_KRR_HCP(CBIG_CODE_DIR, subject_list, unrestricted_csv, restricted_csv, FD_file, FC_file, output_dir)

% CBIG_KRDNN_KRR_HCP(CBIG_CODE_DIR, subject_list, unrestricted_csv, restricted_csv, FD_file, FC_file, output_dir)
% 
% This function runs the kernel ridge regression algorithm to predict
% behavioral measures for HCP dataset. 
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
% 
%   - subject_list
%     Full path of subject list of HCP dataset. It should be a txt file that
%     contains #subject of line, while each line is the subject id of 1
%     subject.
% 
%   - unrestricted_csv
%     Full path of the unrestricted csv for HCP dataset. It should contains the
%     unrestricted behavioral measures.
%
%   - restricted_csv
%     Full path of the restricted csv for HCP dataset. It should contains the
%     restricted behavioral measures.
% 
%   - FD_file
%     Full path of the framewise displacement (FD) for each subject. It should
%     be a txt file that contains #subject of line, while each line
%     corresponds to the mean FD value of a subject. The order should follow
%     the subject order of "subject_list" 
%  
%   - FC_file
%     Full path of the functional connectivity matrix. A matrix "corr_mat" is
%     assumed to be saved in this file. "corr_mat" should be a 3-D matrix 
%     with dimension of #ROI x #ROI x #subjects. Since it is a connectivity
%     matrix, only the lower-triangular off-diagonal entries will be
%     considered as features because the connectivity matrix is symmetric. The
%     order of "corr_mat" should follow the subject order of "subject_list" 
% 
%   - output_dir
%     Full path of the output directory.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

currentdir = pwd;

%% add kernel regression code to path
addpath(genpath(fullfile(CBIG_CODE_DIR, ...
	'/utilities/matlab/predictive_models/KernelRidgeRegression/')));

%% generate setup file
param = load(fullfile(CBIG_CODE_DIR, ...
	'stable_projects/preprocessing/Li2019_GSR/examples/output/KernelRidgeRegression/setup_file.mat'));
% load 20 folds split
temp = load(fullfile(currentdir, 'input', 'kr_hcp', 'He2019_hcp_953_split.mat'));
param.sub_fold = temp.sub_fold;
% output directory
param.outdir = output_dir;
mkdir(param.outdir)
% get matrix of behavioral measures
csv_files = {unrestricted_csv};
subject_header = 'Subject';
measures_set = {'Cognitive','Personality_Task','Social_Emotion'};
y_names = {};
for i = measures_set
	name = i{1};
	measure_list = fullfile(currentdir, 'input', 'kr_hcp', 'measures_lists', [name '_unrestricted.txt']);
	temp = read_sub_list(measure_list);
	y_names = [y_names, temp];
end
y_types = cell(1, size(y_names, 2));
y_types(:) = {'continuous'};
outname = fullfile(output_dir, 'beh_measures.mat');
delimiter = ',';
y = CBIG_read_y_from_csv(csv_files, subject_header, y_names, y_types, ...
	subject_list, outname, delimiter);
param.y = y;
% get matrix of the covariates
csv_files = {unrestricted_csv, restricted_csv};
subject_header = 'Subject';
covariate_names = {'Age_in_Yrs', 'Gender', 'FD'};
covariate_types = {'continuous', 'categorical'};
DVARS_file = 'none';
outname = fullfile(output_dir, 'covariates.mat');
delimiter = ',';
covariates = CBIG_generate_covariates_from_csv(...
	csv_files, subject_header, covariate_names, covariate_types, ...
	subject_list, FD_file, DVARS_file, outname, delimiter);
param.covariates = covariates;
% load FC matrix
temp = load(FC_file);
param.feature_mat = temp.corr_mat;
% number of inner folds
param.num_inner_folds = 20;
% range of lambda for kernel regression
param.lambda_set = [0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 ...
	0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20 30 40 50 60 70 80 100 150 200 ...
	300 500 700 1000 10000 100000 1000000];
save(fullfile(output_dir, 'setup_file.mat'), '-struct', 'param');

%% run kernel regression
CBIG_KRR_workflow(fullfile(output_dir, 'setup_file.mat'), 0);

%% remove kernel regression code to path
rmpath(genpath(fullfile(CBIG_CODE_DIR, ...
	'/utilities/matlab/predictive_models/KernelRidgeRegression/')));

end

function subj_list = read_sub_list(subject_text_list)
% this function will output a 1xN cell where N is the number of
% subjects in the text_list, each subject will be represented by one
% line in the text file
% NOTE: multiple runs of the same subject will still stay on the same
% line
% Each cell will contain the location of the subject, for e.g.
% '<full_path>/subject1_run1_bold.nii.gz <full_path>/subject1_run2_bold.nii.gz'
    fid = fopen(subject_text_list, 'r');
    i = 0;
    while(1);
        tmp = fgetl(fid);
        if(tmp == -1)
            break
        else
            i = i + 1;
            subj_list{i} = tmp;
        end
    end
    fclose(fid);
end