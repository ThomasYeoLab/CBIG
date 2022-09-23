function sub_fold = CBIG_cross_validation_data_split( subject_list, family_csv, ...
    subject_header, family_header, num_folds, seed, outdir, delimiter )

% sub_fold = CBIG_cross_validation_data_split( subject_list, family_csv, ...
%     subject_header, family_header, num_folds, seed, outdir, delimiter )
% 
% This function splits all subjects in "subject_list" into "num_folds"
% folds to perform cross-validation, considering the family structure given
% in "family_csv".
% 
% Inputs:
%   - subject_list
%     Full path of a list of subjects that need to be split. Each line in
%     this list corresponds to one subject ID.
% 
%   - family_csv
%     Full path of the CSV file that contains the family information of all
%     subjects. When all subjects are independent, you can input 'NONE'.
% 
%   - subject_header
%     A string that corresponds to the header of the subject ID column in
%     "family_csv" file. When all subjects are independent, you can input
%     'NONE'.
% 
%   - family_header
%     A string corresponding to the header of the family ID column in
%     "family_csv" file. When all subjects are independent, you can input
%     'NONE'.
% 
%   - num_folds
%     A scalar respresenting how many folds that the subjects need to be
%     split into.
% 
%   - seed
%     A scalar, the random seed used to reorder the families to generate
%     different data splits. In many cases, cross-validation needs to be
%     done multiple times with different data splits to obtain stable
%     results.
% 
%   - outdir
%     Full path of the output directory. This function will output a .mat
%     file called "no_relative_<num_folds>_fold_sub_list.mat" in "outdir".
% 
%   - delimiter
%     The delimiter character used to separate columns in "family_csv".
% 
% Outputs:
%   - sub_fold
%     A struct with two fields: subject_list and fold_index.
%     sub_fold.subject_list is an num_folds x 1 cell, where i-th cell entry
%                           contains the test subject IDs of i-th fold.
%     sub_fold.fold_index is an num_fold x 1 cell, where every entry is a
%                         #subjects x 1 binary vector. 1 means this subject
%                         is assigned to this test fold; 0 means not.
%     This "sub_fold" structure will be saved into "outdir".
% 
% Example:
%   sub_fold = CBIG_cross_validation_data_split( '<path_to_list>/subject_953.txt', ...
%   '<path_to_csv>/RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv', ...
%   'Subject', 'Family_ID', 20, 1, ...
%   '/data/users/jingweil/storage/MyProject/GSR/code_release/test/kernel_regression/HCP', ',' )
% 
% Written by Jingwei, Ru(by) and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(~exist('delimiter', 'var'))
    delimiter = ',';
end

% Read subject IDs and family IDs
[subjects, num_sub] = CBIG_text2cell(subject_list);
fprintf('Input subject list contains %d subjects.\n', num_sub);

if(~strcmpi(family_csv, 'none'))
    assert(~strcmpi(family_header, 'none'), '''family_header'' needs to be passed in.')
    family_IDs = CBIG_parse_delimited_txtfile(family_csv, {family_header}, [], subject_header, subjects, delimiter);
else
    family_IDs = num2cell([1:num_sub]');
    family_IDs = cellfun(@num2str, family_IDs, 'uniformoutput', false);
    warning('No family structure given. Each subject will be treated independently.')
end

% Create family structure
unique_families = unique(family_IDs);
num_fam = length(unique_families);
fprintf('Total number of families: %d.\n', num_fam);

sub_perfam = cell(num_fam, 1);
for i = 1:num_fam
    fam_ind = strcmp(family_IDs, unique_families{i})==1;
    sub_perfam{i} = subjects(fam_ind)';
end

% split families into folds with random seed
fold_list = cell(num_folds,1);
subfold_amount = ceil(num_sub/num_folds);
rng(seed, 'twister');
index = randperm(num_fam);
curr_fold = 0;
for i = 1:num_fam
    curr_fold = mod(curr_fold,num_folds) + 1; 
    while size(fold_list{curr_fold},1)>=subfold_amount
        curr_fold = mod(curr_fold,num_folds) + 1; 
    end
    fold_list{curr_fold} = [fold_list{curr_fold}; sub_perfam{index(i)}];
end

% sort the subject inside fold & create the fold info output
subjects = subjects';
fold_index = cell(num_folds,1);
for i = 1:num_folds
    fold_list{i} = sort(fold_list{i});
    
    fold_index{i} = zeros(num_sub,1);
    for j = 1:size(fold_list{i},1)
        fold_index{i} = fold_index{i} + (strcmp(subjects, fold_list{i}{j}));
    end
    fold_index{i} = logical(fold_index{i});
end
sub_fold = struct('subject_list', fold_list, 'fold_index', fold_index);

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(fullfile(outdir, ['no_relative_' num2str(num_folds) '_fold_sub_list.mat']), 'sub_fold');



end

