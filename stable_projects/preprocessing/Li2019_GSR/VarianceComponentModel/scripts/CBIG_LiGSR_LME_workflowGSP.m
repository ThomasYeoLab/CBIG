function CBIG_LiGSR_LME_workflowGSP( RSFC_file, data_csv, trait_list, covariate_list, ...
    FD_file, DVARS_file, subject_list, outdir, ystem, d, num_samples, rmsub_prefix )

% CBIG_LiGSR_LME_workflowGSP( RSFC_file, data_csv, trait_list, covariate_list, ...
%     FD_file, DVARS_file, subject_list, outdir, ystem, d, num_samples, rmsub_prefix )
% 
% This function performs the whole workflow of the variance component model
% analyses in Li et al., 2019 for a single preprocessing pipeline in the
% Brain Genomics Superstruct Project (GSP) dataset, given a "trait_list".
% Steps include:
% (1) compute FSM from RSFC
% (2) run variance component model on the full subject list
% (3) generate jackknife samples
% (4) run variance component model on each jackknife sample.
% 
% Inputs:
%   - RSFC_file
%     The resting-state functional connectivity filename (a string). It is
%     a .mat file assuming that a #ROIs x #ROIs x #subjects matrix called
%     "corr_mat" is saved in "RSFC_file".
% 
%   - data_csv
%     A string. The full path of the CSV file containing the covariates and
%     traits of all subjects from the GSP dataset.
% 
%   - trait_list
%     A string. The full path of a text file of trait names. Each line
%     in this file corresponds to one trait. The trait names should
%     be consistent with the headers in "data_csv".
%
%   - covariate_list
%     A string. The full path of a text file containing the name of all
%     covariates. Each line corresponds to one covariate. Except for the
%     covariates FD and DVARS, all the other covariate names should exist
%     as headers in "data_csv".
%
%   - FD_file (optional)
%     A string. The full path to a text file of the mean framewise
%     displacement (FD) of all subjects. The number of lines in "FD_file"
%     should be the same as the number of lines in "sub_list". 
%     If the user intends to regress out the covariate 'FD', the user must
%     include 'FD' in the covariate_list. 
%     If the covariates do not include
%     FD, this input variable is not needed. The user can pass in 'NONE'.
% 
%   - DVARS_file (optional)
%     A string. The full path to a text file of the mean DVARS of all
%     subjects. The number of lines in "DVARS_file" should be the same as
%     the number of lines in "subject_list". 
%     If the user intends to regress out the covariate 'DVARS', the user
%     must include 'DVARS' (or 'DV') in the covariate_list.
%     If the covariates do not include DVARS, this input variable is not
%     needed. The user can pass in 'NONE'.
%
%   - subject_list
%     A string. The full path of a text file containing all subject IDs.
%     Each line is one subject ID.
% 
%   - outdir
%     A string. The full path of the output directory.
% 
%   - ystem
%     A string. The trait scores will be read from "data_csv" and saved in
%     a .mat file: [outdir '/y_' ystem '.mat']. For example, if
%     "trait_list" contains 23 behavioral names, you can set ystem =
%     '23behaviors'.
% 
%   - d
%     A scalar or a string, the number of subjects to be removed for each
%     jackknife sample, e.g. 431.
% 
%   - num_samples
%     A scalar or a string, the total number of jackknife samples, e.g.
%     1000.
%   
%   - rmsub_prefix
%     A string, the prefix to the filename of the subject list to be
%     removed for jackknife samples. The list of remove subjects of each
%     jackknife sample will be saved under the output file name
%     [outdir '/jackknife_lists/' rmsub_prefix '_choose' num2str(d) '_set' num2str(i) '.txt']
%     where i ranges from  1 to num_samples.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities'))

if(ischar(d))
    d = str2num(d);
end
if(ischar(num_samples))
    num_samples = str2num(num_samples);
end

%% Compute FSM from RSFC
fprintf('# Step 1. Compute FSM from RSFC.\n')
FSM_file = [outdir '/FSM.mat'];
if(~exist(FSM_file, 'file'))
    CBIG_LiGSR_compute_FSM_from_FC(RSFC_file, FSM_file, 'corr');
end


%% explained trait variance for the full set
fprintf('\n# Step 2. Estimate the explained trait variace on the full set.\n')
CBIG_LiGSR_explained_variance_GSP(data_csv, trait_list, covariate_list, ...
    FD_file, DVARS_file, subject_list, 'none', FSM_file, outdir, ystem, 'fullset');

%% generate jakknife samples
fprintf('\n# Step 3. Generate the list of subjects to be removed for each jackknife sample.\n')
CBIG_LiGSR_NchooseD_families( subject_list, 'none', 'none', 'none', d, num_samples, ...
    fullfile(outdir, 'jackknife_lists'), rmsub_prefix )

%% explained variance for each jakknife sample
fprintf('\n# Step 4. Estimate the explained trait variance on each jackknife sample.\n')
parfor i = 1:num_samples
    rm_sub_list = fullfile(outdir, 'jackknife_lists', [rmsub_prefix '_choose' num2str(d) '_set' num2str(i) '.txt']);
    CBIG_LiGSR_explained_variance_GSP(data_csv, trait_list, covariate_list, ...
        FD_file, DVARS_file, subject_list, rm_sub_list, FSM_file, outdir, ystem, ['del' num2str(d) '_set' num2str(i)], 0);
end

fprintf('\nDone!!!\n')
rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'preprocessing', ...
    'Li2019_GSR', 'VarianceComponentModel', 'scripts', 'utilities'))

end

