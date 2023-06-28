function CBIG_MM_HCP_krr_cv_filter(CBIG_CODE_DIR, HCP_dir, measure_list_dir, output_dir)

% CBIG_MM_HCP_krr_cv_filter(CBIG_CODE_DIR, HCP_dir, measure_list_dir, output_dir)
% 
% This function runs the kernel ridge regression algorithm to each 
% phenotypes (non-brain-imaging phenotypes) for HCP dataset for filter
% phenotypes. 
%
% Inputs:
%   - CBIG_CODE_DIR
%     Full path of the ThomasYeoLab/CBIG repository in your local place
%     (https://github.com/ThomasYeoLab/CBIG).
% 
%   - HCP_dir
%     Path of the phenotypes csv location.
%
%   - measure_list_dir
%     Path of the HCP S1200 phenotypes list directory.
%
%   - output_dir
%     Full path of the output directory for HCP data parse.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add kernel regression code to path
addpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
    'predictive_models', 'KernelRidgeRegression')));

%% generate setup file
param = load(fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing',...
    'Li2019_GSR', 'examples', 'ref_output', 'KernelRidgeRegression',...
    'setup_file.mat'));

output_CV_dir = fullfile(output_dir, 'KRR_CV_output');
if ~exist(output_CV_dir, 'dir')
    mkdir(output_CV_dir);
end

restricted_csv = fullfile(HCP_dir, ['restricted_hc' 'p_data'],...
    'RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv');
subject_list = fullfile(output_dir, 'HCP_diff_roi_subj_list.txt');

% get matrix of behavioral measures
unrestricted_csv = fullfile(HCP_dir, 'subject_measures',...
    'unrestricted_jingweili_12_7_2017_21_0_16_NEO_A_corrected.csv');
csv_files = {unrestricted_csv};
subject_header = 'Subject';
measures_set = {'Cognitive','Personality_Task','Social_Emotion'};
y_names = {};
for i = measures_set
    name = i{1};
    measure_list = fullfile(measure_list_dir, [name '_unrestricted.txt']);
    temp = read_sub_list(measure_list);
    y_names = [y_names, temp];
end
y_types = cell(1, size(y_names, 2));
y_types(:) = {'continuous'};
outname = fullfile(output_CV_dir, 'beh_measures.mat');
delimiter = ',';
y = CBIG_read_y_from_csv(csv_files, subject_header, y_names, y_types, ...
    subject_list, outname, delimiter);
param.y = y;

% only use 1 as covariates
param.covariates = []; % 'none';

% load FC matrix
FC_file = fullfile(output_dir, 'HCP_FC_S1200_1019.mat');
temp = load(FC_file);
param.feature_mat = temp.corr_mat;
% number of inner folds

param.num_inner_folds = 10;
% metric
param.metric = 'none';
% do not save out kernel
param.save_kernel = 0;
% range of lambda for kernel regression
param.lambda_set = [0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 ...
        0.4 0.7 1 1.5 2 2.5 3 3.5 4 5 10 15 20];

rngs = 100;
final_res = zeros(param.num_inner_folds, size(y, 2));

for i = 1:rngs

    % output directory
    output_rng_dir = fullfile(output_CV_dir, ['output_rng_' num2str(i)]);

    res_file = fullfile(output_rng_dir, 'final_result.mat');
    if isfile(res_file)
        disp(['CV with rng ', num2str(i), ' already run, result exists...'])
        res = load(res_file);
        res = res.optimal_acc;
        final_res = final_res + res;
        continue
    end
    param.outdir = output_rng_dir;
    mkdir(param.outdir)

    % 10 folds split
    sub_fold = CBIG_cross_validation_data_split(subject_list, ...
        restricted_csv, 'Subject', 'Family_ID', 10, i, param.outdir, ',' );
    param.sub_fold = sub_fold;

    save(fullfile(output_rng_dir, 'setup_file.mat'), '-struct', 'param');

    %% run kernel regression
    CBIG_KRR_workflow_LITE(fullfile(output_rng_dir, 'setup_file.mat'), 0);
end

final_res = final_res / rngs;
final_beh = mean(final_res)> 0.1;

y_names_filtered = y_names(final_beh);

behavior_list = fullfile(output_dir,...
    'HCP_diff_roi_final_phe_list.txt');
fileID = fopen(behavior_list, 'w');
for i = 1:size(y_names_filtered, 2)
    fprintf(fileID, '%s\n', y_names_filtered{i});
end
fclose(fileID);

y = y(:, final_beh);
phe_corr = final_res;
phe_name = y_names;
phe_corr_filtered = final_res(final_beh);
phe_name_filtered = y_names_filtered;

save(fullfile(output_dir, 'HCP_diff_roi_final_phe.mat'), 'y', 'phe_corr',...
    'phe_name', 'phe_corr_filtered', 'phe_name_filtered');

%% remove kernel regression code to path
rmpath(genpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab',...
    'predictive_models', 'KernelRidgeRegression')));
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
