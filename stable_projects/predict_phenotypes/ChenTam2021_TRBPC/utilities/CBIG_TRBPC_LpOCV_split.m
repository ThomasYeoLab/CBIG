function sub_fold = CBIG_TRBPC_LpOCV_split( subject_list, site_csv, ...
    subject_header, site_header, num_leave_out, outdir, delimiter )

% sub_fold = CBIG_TRBPC_LpOCV_split( subject_list, site_csv, ...
%     subject_header, site_header, num_leave_out, outdir, delimiter )
% 
% This function splits all subjects in "subject_list" into training and
% test folds to perform Leave-p-sites-Out Cross-Validation (LpOCV), each
% split will have subjects from p sites in test set and remaining subjects
% in the training set. The splits will be repeated N-choose-p times. N is
% the total number of sites.
% 
% Inputs:
%   - subject_list
%     Full path of a list of subjects that need to be split. Each line in
%     this list corresponds to one subject ID.
% 
%   - site_csv
%     Full path of the CSV file that contains the site information of all
%     subjects. 
%
%   - subject_header
%     A string that corresponds to the header of the subject ID column in
%     "site_csv" file.
% 
%   - site_header
%     A string corresponding to the header of the site ID column in
%     "site_csv" file.
% 
%   - num_leave_out
%     A scalar. Value of p in the Leave-p-Out Cross-Validation.
% 
%   - outdir
%     Full path of the output directory. This function will output a .mat
%     file called "no_relative_<num_folds>_fold_sub_list.mat" in "outdir".
% 
%   - delimiter
%     The delimiter character used to separate columns in "site_csv".
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
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(~exist('delimiter', 'var'))
    delimiter = ',';
end

% Read subject IDs and site IDs
[subjects, num_sub] = CBIG_text2cell(subject_list);
fprintf('Input subject list contains %d subjects.\n', num_sub);

if(~strcmpi(site_csv, 'none'))
    assert(~strcmpi(site_header, 'none'), '''site_header'' needs to be passed in.')
    site_IDs = CBIG_parse_delimited_txtfile(site_csv, {site_header}, [], subject_header, subjects, delimiter);
else
    site_IDs = num2cell([1:num_sub]');
    site_IDs = cellfun(@num2str, site_IDs, 'uniformoutput', false);
    warning('No site structure given. Each subject will be treated independently.')
end

% Create site structure
unique_sites = unique(site_IDs);
num_site = length(unique_sites);
fprintf('Total number of sites: %d.\n', num_site);

combinations = combnk(unique_sites,num_leave_out);
num_test_fold = size(combinations,1);

fold_list = cell(num_test_fold,1);
fold_index = cell(num_test_fold,1);

for i = 1:num_test_fold
    fold_index{i} = ismember(site_IDs,combinations(i,:));
    fold_list{i} = subjects(fold_index{i});
end

sub_fold = struct('subject_list', fold_list, 'fold_index', fold_index);

if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(fullfile(outdir, ['no_relative_' num2str(num_leave_out) '_fold_sub_list.mat']), 'sub_fold');
fid = fopen(fullfile(outdir,'num_test_folds.txt'),'wt');
fprintf(fid,'%s',num2str(num_test_fold));
fclose(fid);

end
