function [regressors, regressors_CCA] = CBIG_ASDf_genRegressors(sub_info_file, sub_id)
% [regressors, regressors_CCA] = CBIG_ASDf_genRegressors(sub_info_file, sub_id)
% 
% Helper function that generates regressors (age, sex, head motion, sites) for subjects
% indicated by sub_id.
% If sub_id is not specified, the regressors will be for all subjects.
%
% Input:
%     - sub_info_file:
%           .csv file consisting of all subjects' demographics information
%     - sub_id (optional):
%           Nx1 cell array, where N is the number of subjects. IDs of subjects of interest
%
% Output:
%     - regressors:
%           Nx(2+P) matrix, where P is the number of sites associated with
%           the subjects of interest. 1st column is age, 2nd column is sex,
%           3rd column is motion, and 4th to last column are binary indicators
%           for sites.
%     - regressors_CCA:
%           Nx(3+P) matrix. 1st column is age, 2nd column is sex,
%           3rd column is motion, and 4th to last column are binary indicators
%           for sites.
%
% Example:
%     [regressors, regressors_CCA] = CBIG_ASDf_genRegressors('../examples/input/subInfo_est.csv', sub_id)
% 
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% Get sub info
if nargin == 1 || isempty(sub_id)
    [id_sites, ~, id_age, id_sex, id_motion, ~, ~, ~] = CBIG_ASDf_getSubData(sub_info_file);
else
    [id_sites, ~, id_age, id_sex, id_motion, ~, ~, ~] = CBIG_ASDf_getSubData(sub_info_file, sub_id);
end

id = id_sites(:,1);
sites = id_sites(:,2);
sites_unique = unique(sites);

%% Generate regressors
regressors = zeros(length(id),2+length(sites_unique));
regressors_CCA = zeros(length(id),3+length(sites_unique));

regressors(:,1) = cell2mat(id_age(:,2));
regressors(:,2) = cell2mat(id_sex(:,2));
regressors(:,3) = cell2mat(id_motion(:,2));
for j = 1:(length(sites_unique)-1)
    indices = strcmp(sites, sites_unique(j));
    regressors(indices,j+3) = 1;
end

regressors_CCA(:,1) = cell2mat(id_age(:,2));
regressors_CCA(:,2) = cell2mat(id_sex(:,2));
regressors_CCA(:,3) = cell2mat(id_motion(:,2));

for j = 1:length(sites_unique)
    indices = strcmp(sites, sites_unique(j));
    regressors_CCA(indices,j+3) = 1;
end

