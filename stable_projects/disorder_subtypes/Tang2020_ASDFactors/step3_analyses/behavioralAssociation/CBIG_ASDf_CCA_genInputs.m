function [id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores, path_behav_scores, id_all, dx, regressors)
% [id_keep, EB, reg] = CBIG_ASDf_CCA_genInputs(ind_scores, path_behav_scores, id_all, dx, regressors)
% 
% This function gets participants' IDs, exchangeability block (for permutation test
% controlling site) and regressors (i.e., age, sex, motion and sites) for CCA analyses.
%
% Input:
%     - ind_scores:
%           1xP matrix, where P is the number of behavioral scores of interest. 
%           These are the column indicies of the behavioral scores in the .csv file specified by path_behav_scores.
%     - path_behav_scores:
%           Absolute path to the .csv file of all available behavioral scores.
%     - id_all:
%           Nx1 matrix, where N is the number of subjects. This is the IDs of all subjects.
%     - dx:
%           Nx1 matrix. Diagnoses of all subjects. Following ABIDE convention, 1 means ASD, 2 means control.
%     - regressors:
%           NxV matrix, where V is the number of regressors.
% Output:
%     - id_keep:
%           IDs of subjects having all the behavioral scores indicated by ind_scores.
%     - EB:
%           Exchangeability block for generating permutation set using palm_quickperms.m function in PALM package.
%     - reg:
%           Regressors for subjects indicated by id_keep.
%
% Example:
%     [id_keep, EB, reg] = CBIG_ASD_CCA_genInputs(ind_scores,
%     '~/unit_test/data/behavior_scores_654.csv', id_all, dx, regressors)
%     This example retrieves IDs of subjects having the set of behavioral
%     scroes specified by ind_scores. It also generates the exchangeability
%     block EB, and regressors (age, sex, motion and sites). The outputs will be
%     used in CCA analysis.
% 
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variable
if  size(ind_scores,1) > 1
    error('Input argument ''ind_scores'' should be a row vector');
end

if  size(id_all,2) > 1
    error('Input argument ''id_all'' should be a column vector');
end

if  size(dx,2) > 1
    error('Input argument ''dx'' should be a column vector');
end

%% Read behavioral scores
T = readtable(path_behav_scores);
behav_scores = table2cell(T);

%% Get regressors, behavioral scores and IDs of all ASD participants
regressors = regressors(dx == 1,:);
behav_scores = behav_scores(dx == 1,:);
id_asd = id_all(dx == 1);

%% Get IDs of subjects having the set of behavioral scores of interest
sum_nans = sum(cellfun(@isnan,behav_scores(:,ind_scores)),2);
indices = (sum_nans == 0);
id_keep = id_asd(indices);

%% Construct regressors of subjects indicated by id_keep
reg = regressors(indices,1:3); % Age, sex, motion
reg_sites = regressors(indices,4:end); % Sites
ind_extra = (sum(reg_sites,1) == 0);
reg_sites(:,ind_extra) = [];
if length(reg_sites) == 1 % If only left with one site, no need site regressor
    reg_sites = [];
end
reg = [reg reg_sites];

%% Construct the exchangeability block for CCA permutation test
if size(reg,2) > 3 % If more than one site, only allow permutation within sites
    EB = zeros(length(id_keep),3);
    EB(:,1) = -1;
    blk = 1;
    for j = 4:size(reg,2)
        ind = (reg(:,j) == 1);
        EB(ind,2) = blk;
        blk = blk + 1;
    end
    EB(:,3) = (1:length(id_keep))';
else % else, permute freely
    EB = ones(length(id_keep),1);
end


end

