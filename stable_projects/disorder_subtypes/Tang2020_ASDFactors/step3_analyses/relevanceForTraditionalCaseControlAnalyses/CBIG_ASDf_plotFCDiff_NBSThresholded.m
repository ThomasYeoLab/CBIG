function CBIG_ASDf_plotFCDiff_NBSThresholded(sub_info_file, id_factor_file, id_con_file, ADJ,...
    component, corrMat_dir, scale_lim, output_name)
% CBIG_ASDf_plotFCDiff_NBSThresholded(sub_info_file, id_factor_file,
% id_con_file, ADJ, component, corrMat_dir, scale_lim, output_name)
%
% This function plots 419x419 matrix of FC difference between ASD subjects (indicated by id_factor)
% and control subjects (indicated by id_con), thresholded by NBS. If
% output_name is specified, the plot will be saved.
%
% Input:
%     - sub_info_file:
%           .csv file containing all subjects' information, e.g., subjects'
%           IDs, age, sex, sites, mean FD etc.
%     - id_factor_file:
%           .txt file containing ASD subjects' IDs in the subgroup.
%     - id_con_file:
%           .txt file containing control subjects' IDs in the subgroup.
%     - ADJ:
%           A MxM sparse matrix, where M is the number of ROIs. This is the 
%           adjacency matrix from NBS
%     - component:
%           The component to be plotted. In ADJ, nodes belonging to the
%           i-th component should have a value of i.
%     - corrMat_dir:
%           Absolute path to the directory where the subjects' FC matrices
%           are saved
%     - scale_lim:
%           Scale limit to plot the 419x419 matrix. E.g., [-0.2 0.2]
%     - output_name (optional):
%           Full path to the directory to save the plotted figure
%
% Example:
%     CBIG_ASDf_plotFCDiff_NBSThresholded('../examples/input/subInfo_inf.csv',
%     id_factor, id_con, ADJ, 1, corrMat_dir, [-0.1 0.1], '~/example_output/FCdiff')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if size(scale_lim, 1) > 1
    error('Input argument ''scale_lim'' should be a row vector');
end

%% Read subjects' IDs
id_factor = cellstr(num2str());

%% Retrieve all subjects' info
[~, id_dx] = CBIG_ASDf_getSubData(sub_info_file);
id_all = id_dx(:,1);
dx_all = cell2mat(id_dx(:,2));

id_asd_all = id_all(dx_all == 1);
id_con_all = id_all(dx_all == 2);

idx_asd = ismember(id_asd_all, id_factor_file);
idx_con = ismember(id_con_all, id_con_file);

%% Load ASD's FC matrices
load(fullfile(corrMat_dir,'lh2lh_400corrmat_ASD_allSub.mat'));
lh2lh_corr_asd = corr_mat(:,:,idx_asd);
load(fullfile(corrMat_dir,'lh2rh_400corrmat_ASD_allSub.mat'));
lh2rh_corr_asd = corr_mat(:,:,idx_asd);
load(fullfile(corrMat_dir,'rh2rh_400corrmat_ASD_allSub.mat'));
rh2rh_corr_asd = corr_mat(:,:,idx_asd);
load(fullfile(corrMat_dir,'lh2subcor_400corrmat_ASD_allSub.mat'));
lh2subcor_corr_asd = corr_mat(:,:,idx_asd);
load(fullfile(corrMat_dir,'rh2subcor_400corrmat_ASD_allSub.mat'));
rh2subcor_corr_asd = corr_mat(:,:,idx_asd);
load(fullfile(corrMat_dir,'subcor2subcor_400corrmat_ASD_allSub.mat'));
subcor2subcor_corr_asd = corr_mat(:,:,idx_asd);

avgCorr_asd = CBIG_ASDf_indivCorr2avgCorr(lh2lh_corr_asd, lh2rh_corr_asd, rh2rh_corr_asd, ...
    lh2subcor_corr_asd, rh2subcor_corr_asd, subcor2subcor_corr_asd);

%% Load control subjects' FC matrices
load(fullfile(corrMat_dir,'lh2lh_400corrmat_Control_allSub.mat'));
lh2lh_corr_con = corr_mat(:,:,idx_con);
load(fullfile(corrMat_dir,'lh2rh_400corrmat_Control_allSub.mat'));
lh2rh_corr_con = corr_mat(:,:,idx_con);
load(fullfile(corrMat_dir,'rh2rh_400corrmat_Control_allSub.mat'));
rh2rh_corr_con = corr_mat(:,:,idx_con);
load(fullfile(corrMat_dir,'lh2subcor_400corrmat_Control_allSub.mat'));
lh2subcor_corr_con = corr_mat(:,:,idx_con);
load(fullfile(corrMat_dir,'rh2subcor_400corrmat_Control_allSub.mat'));
rh2subcor_corr_con = corr_mat(:,:,idx_con);
load(fullfile(corrMat_dir,'subcor2subcor_400corrmat_Control_allSub.mat'));
subcor2subcor_corr_con = corr_mat(:,:,idx_con);

avgCorr_con = CBIG_ASDf_indivCorr2avgCorr(lh2lh_corr_con, lh2rh_corr_con, rh2rh_corr_con,...
    lh2subcor_corr_con, rh2subcor_corr_con, subcor2subcor_corr_con);

%% Plot FC difference thresholded by NBS
ADJ = full(ADJ);
ADJ(ADJ~=component) = 0; % set other components to 0

corrMat_diff = ADJ.*(avgCorr_asd - avgCorr_con); % ASD - control

if nargin > 8 && ~isempty(output_name)
    CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(corrMat_diff, scale_lim, output_name);
else
    CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(corrMat_diff, scale_lim);
end


end
