function permSet = CBIG_ASDf_genPermSetForNBS(id_factor, id_con, sub_info_file, Nperm)
% permSet = CBIG_ASDf_genPermSetForNBS(id_factor, id_con, sub_info_file, Nperm)
% 
% This function generates the within-site permutation set using PALM package.
% ASD and control subjects from the same site are permuted.
% 
% Input:
%     - id_factor:
%           N1x1 cell array, where N1 is the number of ASD subjects.
%           IDs of ASD subjects.
%     - id_con:
%           N2x1 cell array, where N2 is the number of control subjects.
%           IDs of control subjects.
%     - sub_info_file:
%           Absolute path to the .csv file where all subjects' demographics
%           information is located
%     - Nperm:
%           Number of permutations
% Output:
%      - permSet:
%           (N1+N2)xNperm matrix. Permutation is performed on subjects (ASD
%           + control) from the same site. The 1st column of permSet is the
%           original order.
%
% Example:
%     permSet = CBIG_ASDf_genPermSetForNBS(id_factor, id_con, '../example/input/subInfo_inf.csv', 10000)
% 
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CBIG_CODE_DIR,'external_packages','matlab','non_default_packages','palm','palm-alpha109'));
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

%% Retrieve site information
id_sites_factor = CBIG_ASDf_getSubData(sub_info_file, id_factor);
sites_factor = id_sites_factor(:,2);

id_sites_con = CBIG_ASDf_getSubData(sub_info_file, id_con);
sites_con = id_sites_con(:,2);

sites_names = unique(sites_factor);

%% Construct exchangeability block
nx = length(id_factor);
ny = length(id_con);
EB  = zeros((nx+ny),3); %exchangeability block
EB(:,1) = -1;
sites_all = [sites_factor; sites_con];
for j = 1:length(sites_names)
    indices = strcmp(sites_all, sites_names{j});
    EB(indices,2) = j;
    EB(indices,3) = (1:sum(indices))';
end

%% Generate permutation set using PALM package
permSet = palm_quickperms([],EB,Nperm); % the 1st column is the original order

%% Remove paths
rmpath(fullfile(CBIG_CODE_DIR,'external_packages','matlab','non_default_packages','palm','palm-alpha109'));
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
