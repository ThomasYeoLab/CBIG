function CBIG_TRBPC_feature_transfer_wrapper(outdir)

% CBIG_TRBPC_compute_PFM_average(outdir)
% 
% This function performs the feature transfer prediction in Chen & Tam 2021
% paper. This function assumes certain folder structures. For people
% outside CBIG, they need to modify the paths.
%
% Inputs:
%
%   - outdir
%     output directory for the PFM average files
%
% Outputs:
%   Output prediction accuracy files will be saved in outdir
%
% Written by Jianzhong Chen & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rep_dir = fullfile(getenv('CBIG_REPDATA_DIR'),'stable_projects','predict_phenotypes','ChenTam2022_TRBPC');
feature_file = fullfile(rep_dir,'KRR','feature_files_allFC.mat');
PFM_dir = fullfile(rep_dir,'results', 'PFM', 'KRR', 'allFC');
pred_dir = fullfile(rep_dir,'results', 'KRR', 'allFC');

hypothesis.cog = [9:18,31:36]; % hypothesis-driven cognition
hypothesis.mh = [1,4:6,28,29]; % hypothesis-driven mental health
hypothesis.pers = [19:27]; % hypothesis-driven personality

datadriven.cog = [9:18, 20:21, 31:36]; % data-driven cognition
datadriven.mh = [1, 4:6, 23]; % data-driven mental health
datadriven.pers = [19, 22, 24:29]; % data-driven personality

CBIG_TRBPC_compute_PFM_average(PFM_dir, outdir, hypothesis, datadriven);

domains = {'cog','pers','mh'};
for i = 1:3
    CBIG_TRBPC_feature_transfer_prediction(feature_file, pred_dir, outdir, 'hypothesis', domains{i}, outdir);
    CBIG_TRBPC_feature_transfer_prediction(feature_file, pred_dir, outdir, 'datadriven', domains{i}, outdir);
end

end