function [group_FC, group_FCD_CDF] = ...
    CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, subject_list, session_list, roi_list)

% [group_FC, group_FCD_CDF] = CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(TC_path, subject_list, 
% session_list, roi_list)
% This function generates a group-level FC and FCD from the Alprazolam
% dataset that is used for permutation test. Note that this function requires access to the subject's time
% series data, which are not open to public under the DUA. Thus this
% function is for CBIG internal usage only.
%
% Input:
%   - TC_path: an absolute path to the directory containing subject time series
%   - subject_list: the path to an n-by-1 matrix, each entry is a 5-digit number (n is the number of subjects)
%   - session_list: an n-by-1 matrix, each entry is a binary value. 0 means
%   drug session, 1 means placebo (n is the number of subjects)
%   session
%   - roi_list: a 72-by-1 binary vector, which entry reprensents if an ROI
%   is retained. 1 means the that ROI is retained, 0 means that the ROI is
%   excluded (due to ROI coverage <= 50% or medial wall)
% Output:
%   - group_FC: a 42-by-42 group-level FC used as input to permutation test
%   - group_FCD_CDF: a 10000-by-1 group-level FCD culumative distribution
%   function (CDF) used as input to permutation test
%
% Example:
% [group_FC, group_FCD_CDF] = ...
% CBIG_pFIC_permutation_test_Alprazolam_group_level_FC_FCD(['isilon/CSC2/Yeolab/Data/Alpraz/TASK_fmriprep/' ... 
%   'desikan_from_vol/'], 'replication/Alprazolam/input/drug/training_subject_list.txt', randi([0, 1], [1, 15])', ...
%   'replication/Alprazolam/input/drug/roi_list_0.5.txt')
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

subject_list = dlmread(subject_list);
%% exclude medial wall
fid = fopen(roi_list);
roi = textscan(fid, '%f');
fclose(fid);
roi = roi{1,1};
roi([1, 5, 37, 41]) = [];
num_roi = sum(roi, 1);

%% initialize group-level FC and FCD_CDF
FC_all = zeros(num_roi, num_roi, length(subject_list));
FCD_all = zeros(length(subject_list), 10000);

% because under the permutation test setup, the session label (placebo vs. drug) is randomly
% shuffled within a subject, thus, the session_list is a random list of 0s
% and 1s
for i = 1:length(subject_list)
    subject = subject_list(i);
    [FC, FCD_CDF] = CBIG_pFIC_generate_Alprazolam_subject_level_FC_FCD(TC_path, subject, session_list(i), roi_list);
    FC_all(:, :, i) = FC;
    FCD_all(i, :) = FCD_CDF;
end

% group-level FC requires Fisher transformation before averaging
FC_sum = zeros(num_roi);
for i = 1:size(FC_all, 3)
   FC = FC_all(:, :, i); 
   FC_sum = FC_sum + CBIG_StableAtanh(FC);   
end
group_FC = tanh(FC_sum/size(FC_all, 3));

group_FCD_CDF = mean(FCD_all);
group_FCD_CDF = round(group_FCD_CDF');

end
