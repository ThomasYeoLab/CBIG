function CBIG_TRBPC_replicate_matrix_plots(outdir)
% CBIG_TRBPC_replicate_matrix_plots(outdir)
% 
% This function will generate all matrix-related figures for the TRBPC project
%
% Inputs:
%
%   - outdir
%     output directory for the figures
%
% Outputs:
%   Output figures and .mat files will be stored in outdir
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

data_dir = fullfile(getenv('CBIG_REPDATA_DIR'),'stable_projects',...
    'predict_phenotypes','ChenTam2022_TRBPC','results');
significant_edges = fullfile(data_dir, 'PFM','KRR','allFC', 'PFM_significance_masks.mat');

datadriven_ind.cog = [9:18, 20:21, 31:36]; % data-driven cognition
datadriven_ind.mh = [1, 4:6, 23]; % data-driven mental health
datadriven_ind.pers = [19, 22, 24:29]; % data-driven personality

hypothesis_ind.cog = [9:18,31:36]; % hypothesis-driven cognition
hypothesis_ind.mh = [1,4:6,28,29]; % hypothesis-driven mental health
hypothesis_ind.pers = [19:27]; % hypothesis-driven personality

CBIG_TRBPC_matrix_plots_wrapper(data_dir, outdir, hypothesis_ind, datadriven_ind, significant_edges);