function [id_sites, id_dx, id_age, id_sex, id_motion, id_fiq, id_scores, id_factorLoading] = CBIG_ASDf_getSubData(...
    sub_info_file,id_sub, behav_score_file, idx_scores, factorLoading_file)
% [id_sites, id_dx, id_age, id_sex, id_motion, id_fiq, id_scores,
% id_factorLoading] = CBIG_ASDf_getSubData(sub_info_file, id_sub, behav_score_file, idx_scores , factorLoading_file)
%
% Helper function that extracts subjects' information from input .csv files. It
% assumes that the columns in sub_info_file are as follows:
%       1st column: Sites
%       2nd column: Subjects' IDs
%       3rd column: Diagnoses
%       4th column: Age
%       5th column: Sex
%       6th column: Head motion
%       7th column: FIQ scores
%
% If only sub_info_file is specified, it will retrieve information of all subjects;
% If id_sub is specified, it will retrieve information of subjects of interest indicated by id_sub;
% If behav_score_file & idx_scores are specified, it will also retrieve behavioral scores;
% If factorLoading_file is entered, it will also retrieve the factor compositions.
%
% Input:
%     - sub_info_file:
%           .csv file containing all subjects' demographics information
%     - id_sub (optional):
%           Nx1 cell array, where N is the number of subjects. IDs of
%           subjects of interest. The ordering of subjects should follow
%           the ordering in sub_info_file.
%     - behav_score_file (optional):
%           .csv file containing all subjects' behavioral data. The odering
%           of subjects should follow the ordering in sub_info_file.
%     - idx_scores (optional):
%           1xQ vector. Indices of behavioral scores of interest. These are
%           the column indices in behav_score_file.
%     - factorLoading_file (optional):
%           .txt file containing factor compositions of ASD subjects
% Output:
%     - id_dx:
%           Nx2 cell array, where N is the number of subjects. 1st column is
%           subjects' IDs, and 2nd column is diagnoses of the subjects.
%     - id_sites:
%           Nx2 cell array, where N is the number of subjects. 1st column is
%           subjects' IDs, and 2nd column is acquisition sites of the subjects.
%     - id_age:
%           Nx2 cell array, where N is the number of subjects. 1st column is
%           subjects' IDs, and 2nd column is age of the subjects.
%     - id_sex:
%           Nx2 cell array, where N is the number of subjects. 1st column is
%           subjects' IDs, and 2nd column is sex of the subjects.
%     - id_motion:
%           Nx2 cell array, where N is the number of subjects. 1st column is
%           subjects' IDs, and 2nd column is mean FD of the subjects.
%     - id_fiq:
%           Nx2 cell array, where N is the number of subjects. 1st column is
%           subjects' IDs, and 2nd column is FIQ score of the subjects.
%     - id_scores:
%           Nx(Q+1) cell array, where Q is the number of behavioral scores of
%           interest. 1st column is subjects' IDs, and 2nd to last columns
%           are behavioral scroes. [] if idx_scores or behav_score_file is not specified.
%     - id_factorLoading:
%           NxK cell array, where K is the number of factors. Factor
%           compositions of subjects of interest. 1st column is subjects'
%           IDs, 2nd to last columns are factor compositions. [] if
%           factorLoading_file is not specified.
%
% Example:
%       [id_sites, id_dx] = CBIG_ASDf_getSubData('../examples/input/subInfo_est.csv')
%       This function will retrieve the IDs, sites & diagnoses of all subjects
%       listed in the file '../examples/input/subInfo_est.csv'.
%
%       [~, id_dx, ~, ~, ~, ~, id_scores] =
%       CBIG_ASDf_getSubData('../examples/input/subInfo_est.csv', id_sub,
%       '../examples/input/behavior_scores.csv', [5 6])
%       This funtion will retrieve the IDs, diagnoses as well as two behavioral
%       scores of subjects specified by id_sub.
%
%       [~, ~, ~, ~, ~, ~, ~, id_factorLoading] =
%       CBIG_ASDf_getSubData('../examples/input/subInfo_est.csv', [], [],
%       [], '~/example_output/visualizeFactors/k3/r94/factorComp.txt')
%       This function will retreive the IDs & factor compositions of all subjects.
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if exist('idx_scores', 'var')
    if size(idx_scores,1) > 1
        error('Input argument ''idx_scores'' should be a row vector');
    end
end

subInfo = readtable(sub_info_file);

id_all = table2cell(subInfo(:,2));
dx_all = table2array(subInfo(:,3));

if ~ischar(id_all{1})
    id_all = cellfun(@num2str,id_all,'UniformOutput',false);
end

if (nargin >= 2) && ~isempty(id_sub)
    id_sub = cellfun(@num2str,id_sub,'UniformOutput',false);
    idx_sub = ismember(id_all, id_sub);
else
    idx_sub = ismember(id_all, id_all);
    id_sub = id_all;
end


% ID and sites
id_sites = [id_sub table2cell(subInfo(idx_sub, 1))]; %column 2 is ID, column 1 is sites

% ID and diagnosis
id_dx = [id_sub table2cell(subInfo(idx_sub, 3))]; %column 2 is ID, column 3 is diagnosis

% ID and age
id_age = [id_sub table2cell(subInfo(idx_sub, 4))]; %column 2 is ID, column 4 is age

% ID and sex
id_sex = [id_sub table2cell(subInfo(idx_sub, 5))]; %column 2 is ID, column 5 is sex

% ID and FD
id_motion = [id_sub table2cell(subInfo(idx_sub, 6))]; %column 2 is ID, column 6 is motion

% ID and FIQ
id_fiq = [id_sub table2cell(subInfo(idx_sub, 7))]; %column 2 is ID, column 7 is FIQ

% ID and behavioral scores
if (nargin >= 4) && ~isempty(idx_scores) && ~isempty(behav_score_file)
    behavScore = readtable(behav_score_file);
    id_scores = [id_sub table2cell(behavScore(idx_sub, idx_scores))];
else
    id_scores = [];
end

%% If the user wants to retrieve factor compositions as well
if (nargin == 5) && (~isempty(factorLoading_file))
    idx_asd_fl = ismember(id_all(dx_all==1,1), id_sub); % Dimension is number of ASD subjects
    idx_fl = ismember(id_all, id_sub); % Dimension is number of all subjects
    factorComp = dlmread(factorLoading_file);
    id_factorLoading = [id_all(idx_fl&(dx_all==1),1) num2cell(factorComp(idx_asd_fl,:))];
else
    id_factorLoading = [];
end
