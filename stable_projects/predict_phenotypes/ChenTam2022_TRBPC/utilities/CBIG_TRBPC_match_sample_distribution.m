function [id_sample, p_before, p_after] = CBIG_TRBPC_match_sample_distribution(...
    cov_table,id_all, id_subset, cov_to_match, cov_types, cov_levels, target_sample)

% [id_sample, p_before, p_after] = CBIG_TRBPC_match_sample_distribution(...
%    cov_table,id_all, id_subset, cov_to_match, cov_types, cov_levels, target_sample)
%
% This fucntion matches the covariates distribution between included and
% excluded subjects using stratified subsampling.
%
% Inputs:
%   - cov_table
%     A table that contains all covariates of all subjects.
%
%   - id_all
%     A cell array contains the subject id of the population we want to
%     match the covariates distribution.
%
%   - id_subset
%     A subset of id all. This function would select subjects from
%     id_subset so that covariates distribution of the selected sample
%     matches the distribution in id_all.
%
%   - cov_to_match
%     A 1*D cell array contains the name of variables that we want to match.
%     D is the number of covariates. The names must be found in the column
%     names in cov_table.
%
%   - cov_types
%     A 1*D cell array while D is the number of covariates. Each cell should be
%     either 'continous' or 'categorical', indicating the variable type.
%
%   - cov_levels
%     A 1*D cell array while D is the number of covariates. Each cell contains the
%     bins we want divide the variable into. For example, for age the
%     levels could be a {[10,20],[20,30]}, so we will divide subjects into 2 bins.
%     age 10 ~ 20 and age 20 ~ 30. For a categorical variable like sex,
%     it should be {'M','F'}, so we divide subjects into 2 bins (male vs
%     female).
%
%   - target_sample
%     An integer indicating the target number of subjects we want to match.
%
% Outputs:
%   - id_sample
%     Samples selected from id_subset that has the same covariates
%     distribution as id_all
%
%   - p_before
%     A 1*D array contains the p-value of covariates distribution between
%     included subjects (id_subset) and excluded subjects
%     (id_all-id_subset). For continous variables, we use t-test. For
%     categorical variables, we use chi-square test.
%
%   - p_after
%     A 1*D array contains the p-value of covariates distribution between
%     included subjects (id_sample) and excluded subjects
%     (id_all-id_sample).
%
% Example: [id_sample, p_before, p_after] = CBIG_TRBPC_match_sample_distribution(cov_table, ...
%    id_all, id_subset, {'age','sex'}, {'continous','categorical'}, {{[10,20],[20,30]},{'M','F'}}, 500)
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rng(1)
N_var = length(cov_to_match);
all_combs = allcomb(cov_levels);
ind = CBIG_find_cell_in_cell(id_subset, id_all);
subset_table = cov_table(ind,:);
id_sample = [];
for i = 1:length(all_combs)
    all_pass_flag = ones(length(id_all),1);
    subset_pass_flag = ones(length(id_subset),1);
    for j = 1:N_var
        curr_var_level = all_combs{i,j};
        curr_var = cov_to_match{j};
        var_type = cov_types{j};
        if strcmp(var_type,'continous')
            curr_all_pass = and(cov_table.(curr_var) >= curr_var_level(1),cov_table.(curr_var) < curr_var_level(2));
            curr_subset_pass = and(subset_table.(curr_var) >= curr_var_level(1),...
            subset_table.(curr_var) < curr_var_level(2));
            
        else
            curr_all_pass = strcmp(cov_table.(curr_var),curr_var_level);
            curr_subset_pass = strcmp(subset_table.(curr_var),curr_var_level);
        end
        all_pass_flag = and(all_pass_flag,curr_all_pass);
        subset_pass_flag = and(subset_pass_flag, curr_subset_pass);
    end
    all_prop = mean(all_pass_flag);
    target_subset_sample = ceil(all_prop*target_sample);
    
    if target_subset_sample < sum(subset_pass_flag)
        id_sample = [id_sample; datasample(id_subset(subset_pass_flag), target_subset_sample, 'Replace', false)];
    else
        id_sample = [id_sample; id_subset(subset_pass_flag)];
    end
end

p_before = zeros(N_var,1);
p_after = zeros(N_var,1);
ind_included_before = contains(id_all,id_subset);
ind_included_after = contains(id_all,id_sample);
for i = 1:N_var
    curr_var = cov_to_match{i};
    if strcmp(cov_types{i},'continous')
        [~,p_before(i)] = ttest2(cov_table.(curr_var)(ind_included_before),cov_table.(curr_var)(~ind_included_before));
        [~,p_after(i)] = ttest2(cov_table.(curr_var)(ind_included_after),cov_table.(curr_var)(~ind_included_after));
    else
        [~,~,p_before(i)] = crosstab(cov_table.(curr_var),ind_included_before);
        [~,~,p_after(i)] = crosstab(cov_table.(curr_var),ind_included_after);
    end
end