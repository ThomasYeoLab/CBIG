function CBIG_KRDNN_KRR_UKBB(setup_file, num_val_tes, lambda_set)

% CBIG_KRDNN_KRR_UKBB(setup_file, num_val_tes, lambda_set)
% 
% This function runs the kernel ridge regression algorithm to predict
% behavioral measures for UK Biobank dataset. 
%
% Inputs:
%   - setup_file
%     Path of the setup file (.txt file).
%     setup_file has the following fields:
%     1. Line 1
%        `data`, the location of your output and intermediate values. You can
%        also change this to any place you want.
%     2. Line 2
%        Path of subject list of UK biobank dataset. It should be a txt file
%        that contains #subject of line (8868 in our paper), while each line
%        is the subject id of 1 subject.
%     3. Line 3
%        Path of CSV file extracted from UK biobank dataset. You can
%        remove the unused measures from the csv to increase the loading
%        speed.
%     4. Line 4
%        Path of the mean rfMRI head motion (data-field 25741) for each
%        subject. It should be a txt file that contains #subject of line,
%        while each line corresponds to the mean rfMRI head motion of a
%        subject. The order should follow the subject order of "subject_list"
%     5. Line 5
%        Path of the functional connectivity matrix. A matrix "corr_mat"
%        is assumed to be saved in this file. "corr_mat" should be a 3-D
%        matrix with dimension of #ROI x #ROI x #subjects. Since it is a
%        connectivity matrix, only the lower-triangular off-diagonal entries
%        will be considered as features because the connectivity matrix is
%        symmetric. The order of "corr_mat" should follow the subject order
%        of "subject_list" 
%     6. Line 6
%        Path of measure list text file
%     7. Line 7
%        Name of the measures lists text file in measure_dir folder. You
%        do not need to change this unless you want to change the measure list
%   - num_val_tes
%     Number (integer) of subjects for validation and testing set.
%   - lambda_set
%     1xN vector, where N is number of lambda to search. This is the set of
%     lambda for kernel regression algorithm to search through.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    

% process input
if ~isnumeric(num_val_tes)
    num_val_tes = str2double(num_val_tes);
end
if ~isscalar(num_val_tes)
    num_val_tes = num_val_tes(1);
end
if floor(num_val_tes) ~= num_val_tes
    num_val_tes = floor(num_val_tes);
end
assert(isvector(lambda_set), 'lambda_set should be a vector')

%% data prepare
% split data into training\validation\test set
CBIG_KRDNN_split_subject(setup_file, 140, num_val_tes)

% regress out age, gender and motion from behaviour measures for each fold
CBIG_KRDNN_regress_from_behaviour(setup_file)

% generate kernel
CBIG_KRDNN_generate_K_train(setup_file)
CBIG_KRDNN_generate_K_valid(setup_file)
CBIG_KRDNN_generate_K_test(setup_file)

%% actual running
% Training
CBIG_KRDNN_train_valid_krr(setup_file, lambda_set)

% Testing
CBIG_KRDNN_test_krr(setup_file, lambda_set)

% Get final result
CBIG_KRDNN_get_ukbb_krr_result(setup_file, lambda_set);
