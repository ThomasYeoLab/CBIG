function permutation_test_input_struct = ...
    CBIG_pFIC_Alprazolam_permutation_test_wrapper(TC_path, training_subject_list, validation_subject_list, ...
    test_subject_list, roi_list)

% permutation_test_input_struct = ...
%   CBIG_pFIC_Alprazolam_permutation_test_wrapper(TC_path, training_subject_list, validation_subject_list, ...
%   test_subject_list, roi_list)
% This function generates 100 sets of FCs and FCD CDFs for permutation test. Note that this function requires 
% access to the subject's time series data, which are not open to public under the DUA. Thus this
% function is for CBIG internal usage only.
%
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - training_subject_list: the path to an 15-by-1 matrix, each entry is a
%   5-digit subject id from the training set
%   - validation_subject_list: the path to an 15-by-1 matrix, each entry is a
%   5-digit subject id from the validation set
%   - test_subject_list: the path to an 15-by-1 matrix, each entry is a
%   5-digit subject id from the test set
%   - roi_list: a 72-by-1 binary vector, which entry reprensents if an ROI
%   is retained. 1 means the that ROI is retained, 0 means that the ROI is
%   excluded (due to ROI coverage <= 50% or medial wall)
% Output:
%   - permutation_test_input_struct: a structure with 12 fields, 6 fields
%   of FCs and 6 fields of FCD CDFs. Both FCs and FCD CDFs are further
%   divided into training, validation and test set of 'drug' and 'placebo'
%   sessions. The index of the last dimension of FC and FCD CDF (3rd
%   dimension for FC, 2nd dimension for FCD) establish the correspondance between the 2.
%
% Example:
% permutation_test_input_struct = ...
%   CBIG_pFIC_Alprazolam_permutation_test_wrapper(['isilon/CSC2/Yeolab/Data/Alpraz/TASK_fmriprep/' ... 
%   'desikan_from_vol/'], ...
%   'replication/Alprazolam/input/drug/training_subject_list.txt', ... 
%   'replication/Alprazolam/input/drug/validation_subject_list.txt', ...
%   'replication/Alprazolam/input/drug/test_subject_list.txt', ...
%   'replication/Alprazolam/input/drug/roi_list_0.5.txt');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% initialization
FC_train_drug = zeros(42, 42, 100);
FC_validation_drug = zeros(42, 42, 100);
FC_test_drug = zeros(42, 42, 100);
FC_train_placebo = zeros(42, 42, 100);
FC_validation_placebo = zeros(42, 42, 100);
FC_test_placebo = zeros(42, 42, 100);

FCD_CDF_train_drug = zeros(10000, 100);
FCD_CDF_validation_drug = zeros(10000, 100);
FCD_CDF_test_drug = zeros(10000, 100);
FCD_CDF_train_placebo = zeros(10000, 100);
FCD_CDF_validation_placebo = zeros(10000, 100);
FCD_CDF_test_placebo = zeros(10000, 100);

%% generate 100 sets permutation test input FCs and FCDs
for seed = 1:100
    disp(['Preparing input FC and FCD CDF for permutation test [' num2str(seed) '] ...'])
    rng(10000+seed, 'twister');
    % first column represents the 'drug' sesion after permutation
    % second cloum represents the 'placebo' sesion after permutation
    % 0 -> drug, 1 -> placebo
    random_session_list = zeros(45, 2);
    random_session_list(:, 1) = transpose(randi([0, 1], [1, 45]));
    random_session_list(:, 2) = 1 - random_session_list(:, 1);   
%% drug
    [FC_train_drug(:, :, seed), FCD_CDF_train_drug(:, seed)] = ...
        CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, training_subject_list, ...
        random_session_list(1:15, 1), roi_list);
    [FC_validation_drug(:, :, seed), FCD_CDF_validation_drug(:, seed)] = ...
        CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, validation_subject_list, ...
        random_session_list(16:30, 1), roi_list);
    [FC_test_drug(:, :, seed), FCD_CDF_test_drug(:, seed)] = ...
        CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, test_subject_list, ...
        random_session_list(31:45, 1), roi_list); 
%% placebo
    [FC_train_placebo(:, :, seed), FCD_CDF_train_placebo(:, seed)] = ...
        CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, training_subject_list, ...
        random_session_list(1:15, 2), roi_list);
    [FC_validation_placebo(:, :, seed), FCD_CDF_validation_placebo(:, seed)] = ...
        CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, validation_subject_list, ...
        random_session_list(16:30, 2), roi_list);
    [FC_test_placebo(:, :, seed), FCD_CDF_test_placebo(:, seed)] = ...
        CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, test_subject_list, ...
        random_session_list(31:45, 2), roi_list);
end

%% prepare output
permutation_test_input_struct.FC_train_drug = FC_train_drug;
permutation_test_input_struct.FC_validation_drug = FC_validation_drug;
permutation_test_input_struct.FC_test_drug = FC_test_drug;
permutation_test_input_struct.FC_train_placebo = FC_train_placebo;
permutation_test_input_struct.FC_validation_placebo = FC_validation_placebo;
permutation_test_input_struct.FC_test_placebo = FC_test_placebo;

permutation_test_input_struct.FCD_CDF_train_drug = FCD_CDF_train_drug;
permutation_test_input_struct.FCD_CDF_validation_drug = FCD_CDF_validation_drug;
permutation_test_input_struct.FCD_CDF_test_drug = FCD_CDF_test_drug;
permutation_test_input_struct.FCD_CDF_train_placebo = FCD_CDF_train_placebo;
permutation_test_input_struct.FCD_CDF_validation_placebo = FCD_CDF_validation_placebo;
permutation_test_input_struct.FCD_CDF_test_placebo = FCD_CDF_test_placebo;
end
