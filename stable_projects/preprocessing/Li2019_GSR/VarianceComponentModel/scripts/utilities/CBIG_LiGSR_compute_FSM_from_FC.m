function CBIG_LiGSR_compute_FSM_from_FC(RSFC_file, FSM_file, type, scale)

% CBIG_LiGSR_compute_FSM_from_FC(RSFC_file, FSM_file, type, scale)
%
% Given the ROIs to ROIs functional connectivity file, this function first
% grab the lower triangle part and reshape the matrix then calls
% CBIG_Kernel_with_scale_factor.m to compute Functional Similarity Matrix
% (FSM). This FSM has the dimension of #subjects x #subjects, and will be
% saved as a .mat file with a name provided by RSFC_file input argument.
%
% Input:
%   - RSFC_file:
%     A string. The full path of a .mat file. This file contains a #ROIs x
%     #ROIs x #subjects matrix called "corr_mat". Entry (i,j,k) in this
%     matrix is the functional connectivity between ROI i and ROI j of
%     subject k.
% 
%   - FSM_file:
%     A string. The output file name of FSM matrix (full path). A matrix
%     named FSM will be saved in this file with dimention of #subjects x
%     #subjects.
%
%   - type:
%     A string. The method to compute FSM. Choose from 'corr',
%     'Exponential', or 'Gaussian'.
%
%   - scale:
%     A scalar. The scaling factor of kernel. If 'corr' is chosen, the user
%     does not need to pass in this parameter.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


load(RSFC_file);
corr_size = size(corr_mat);

% lower triangluar entries, reshape
tril_ind = find(tril(ones(corr_size(1:2)), -1));
corr_mat_reshape = reshape(corr_mat, corr_size(1)*corr_size(2), corr_size(3));
feature = corr_mat_reshape(tril_ind, :);

% compute FSM
if(strcmp(type, 'corr'))
    FSM = CBIG_crossvalid_kernel_with_scale(feature, [], [], [], type);
else
    FSM = CBIG_crossvalid_kernel_with_scale(feature, [], [], [], type, scale);
end

[FSM_path, ~, ~] = fileparts(FSM_file);
mkdir(FSM_path);
save(FSM_file, 'FSM', '-v7.3');



end