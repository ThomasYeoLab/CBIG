function CBIG_MM_KRR_filter_split_subject(data_dir, subject_list_file, rng_num, per_valid_test)

% CBIG_MM_KRR_filter_split_subject(data_dir, subject_list_file, rng_num, per_valid_test)
% 
% This function splits subjects into training, validation and testing set.
%
% Inputs:
%   - data_dir
%     Path of the your output and intermediate values. You can
%     also change this to any place you want.
%   - subject_list_file
%     Path of subject list of setup_file dataset. It should be a txt file
%     that contains #subject of line, while each line is the subject id of
%   - rng_num
%     Number (integer) of random number generator
%   - num_val_tes
%     percentage (float) of subjects for validation and testing set.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rng(rng_num);

% read out subject file
fileID = fopen(subject_list_file);
temp = textscan(fileID, '%d');
subject_list = temp{1};

% Train/Validation/Test split
num_total = size(subject_list, 1);
num_valid = floor(num_total * per_valid_test);
num_test = floor(num_total * per_valid_test);
num_train = num_total - num_test - num_valid;

% generate list
subject_list_rand = subject_list(randperm(num_total));
subject_list_valid = sort(subject_list_rand(1:num_valid));
subject_list_test = sort(subject_list_rand(num_valid+1:num_valid+num_test));
subject_list_train = sort(subject_list_rand(end-num_train+1:end));

% generate index list
fold_index = zeros(num_total,1);
for i = 1:size(subject_list_valid,1)
    fold_index = fold_index + (subject_list_valid(i) == subject_list);
end
fold_index = fold_index * 2;
for i = 1:size(subject_list_test,1)
    fold_index = fold_index + (subject_list_test(i) == subject_list);
end

% save them in strucuture
sub_fold = struct('fold_index', fold_index);

% save out
save(fullfile(data_dir, 'ukbb_subject_split.mat'), 'sub_fold');

end