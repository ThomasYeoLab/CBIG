function CBIG_MM_KRR_MM_split_subject(data_dir, subject_list_file, rng_num, per_valid,...
    subj_meta_train, subj_meta_test, pre_fix)

% CBIG_MM_KRR_MM_split_subject(data_dir, subject_list_file, rng_num, per_valid, subj_meta_train, subj_meta_test)
% 
% This function splits subjects into k subjects and remaining subjects.
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
%   - per_valid
%     Portion or percentage (float) of validation data in total data
%   - subj_meta_train
%     Full path of subject list of training meta-set. It should be a txt
%     file that contains #subject of line, while each line is the subject
%     id of 1 subject.
%   - subj_meta_test
%     Full path of subject list of test meta-set. It should be a txt
%     file that contains #subject of line, while each line is the subject
%     id of 1 subject.
%   - pre_fix
%     name pre fix for the split file saving
%
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rng(rng_num);

% read out subject file
fileID = fopen(subject_list_file);
temp = textscan(fileID, '%d');
subject_list = temp{1};

fileID = fopen(subj_meta_train);
temp = textscan(fileID, '%d');
subj_meta_train = temp{1};

fileID = fopen(subj_meta_test);
temp = textscan(fileID, '%d');
subj_meta_test = temp{1};

% Train/Validation/Test split
num_total = size(subj_meta_train, 1);
num_valid = floor(num_total * per_valid);
num_train = num_total - num_valid;

% generate list
subject_list_rand = subj_meta_train(randperm(num_total));
subject_list_train = sort(subject_list_rand(1:num_train));
subject_list_valid = sort(subject_list_rand(end-num_valid+1:end));
subject_list_test = subj_meta_test;

tmp = sort(cat(1, subject_list_valid, subject_list_train, subject_list_test));
if ~isequal(tmp, subject_list)
    disp(tmp)
    error('subject not match during split')
end

% generate index list
fold_index = zeros(size(subject_list, 1), 1);
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
save(fullfile(data_dir, [pre_fix '_subject_split.mat']), 'sub_fold');

end