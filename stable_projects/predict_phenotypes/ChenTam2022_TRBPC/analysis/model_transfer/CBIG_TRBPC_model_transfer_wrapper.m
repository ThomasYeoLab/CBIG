function CBIG_TRBPC_model_transfer_wrapper(outdir)

% CBIG_TRBPC_model_transfer_wrapper(outdir)
% 
% This function performs the model transfer prediction in Chen & Tam 2021
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
pred_dir = fullfile(rep_dir,'results', 'KRR', 'allFC');

hypothesis_ind.cog = [9:18,31:36]; % hypothesis-driven cognition
hypothesis_ind.mh = [1,4:6,28,29]; % hypothesis-driven mental health
hypothesis_ind.pers = [19:27]; % hypothesis-driven personality

datadriven_ind.cog = [9:18, 20:21, 31:36]; % data-driven cognition
datadriven_ind.mh = [1, 4:6, 23]; % data-driven mental health
datadriven_ind.pers = [19, 22, 24:29]; % data-driven personality

CBIG_TRBPC_model_transfer(pred_dir, hypothesis_ind, datadriven_ind, outdir)