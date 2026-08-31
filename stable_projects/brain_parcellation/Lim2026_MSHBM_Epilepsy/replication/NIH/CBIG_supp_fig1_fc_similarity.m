% ====================================
%% Help
% ====================================
%
% CBIG_supp_fig1_fc_similarity
%
% This script replicates Supplementary Figure 1: the vertex-wise Pearson
% correlation between the NIH-34 and HCP-40 group-average whole-brain
% functional connectivity (FC) fingerprint maps. For each cortical vertex,
% the FC profile (its row of the group-average FC matrix) is correlated
% between the two cohorts, producing a similarity map displayed on the
% fsaverage6 inflated surface.
%
% PREREQUISITES:
%   Pre-computed group-average FC block files for NIH-34 and HCP-40 must
%   exist in the input directories below. Each block file contains:
%     group_FC_avg: (n_vertices x block_size) FC matrix block
%     row_idx:      vertex indices corresponding to the block columns
%
% INPUTS:
%   CBIG_CODE_DIR: environment variable pointing to the CBIG codebase root.
%   NIH-34 FC blocks: replication/NIH/input/FC_compare/NIH_34_blocks/
%     FC_block_001.mat ... FC_block_082.mat
%   HCP-40 FC blocks: replication/NIH/input/FC_compare/HCP_group/
%     FC_block_001.mat ... FC_block_082.mat
%
% OUTPUTS:
%   Prints mean similarity across cortical vertices to console.
%   Displays Supplementary Figure 1 on screen (not saved).
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Last tested on 20 May 2026.

% ====================================
%% Setup Paths
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
proj_dir = fullfile(CBIG_CODE_DIR, ...
    'stable_projects/brain_parcellation/Lim2026_MSHBM_Epilepsy');
data_dir = getenv('CBIG_EPILEPSY_DATA_DIR');

addpath(fullfile(CBIG_CODE_DIR, 'utilities/matlab/figure_utilities'));

nih_block_dir = fullfile(data_dir, 'NIH/input/FC_compare/NIH_34_blocks');
hcp_block_dir = fullfile(data_dir, 'NIH/input/FC_compare/HCP_group');

% ====================================
%% Compute FC Similarity Map
% ====================================

n_vertices_total = 81924;
block_size       = 1000;
n_blocks         = ceil(n_vertices_total / block_size);

FC_similarity_map = NaN(n_vertices_total, 1);

fprintf('\n====== Computing NIH-HCP FC similarity map ======\n');

for block_idx = 1:n_blocks
    fprintf('  Block %d / %d\n', block_idx, n_blocks);

    d_NIH = load(fullfile(nih_block_dir, sprintf('FC_block_%03d.mat', block_idx)));
    d_HCP = load(fullfile(hcp_block_dir, sprintf('FC_block_%03d.mat', block_idx)));

    FC_NIH  = d_NIH.group_FC_avg;
    FC_HCP  = d_HCP.group_FC_avg;
    row_idx = d_NIH.row_idx;

    for c = 1:length(row_idx)
        v           = row_idx(c);
        nih_col     = FC_NIH(:, c);
        hcp_col     = FC_HCP(:, c);
        valid_rows  = ~isnan(nih_col) & ~isnan(hcp_col);
        if sum(valid_rows) < 2; continue; end
        FC_similarity_map(v) = corr(nih_col(valid_rows), hcp_col(valid_rows));
    end
end

fprintf('Mean similarity (cortical vertices): %.4f\n', ...
    mean(FC_similarity_map(~isnan(FC_similarity_map))));

% ====================================
%% Plot Supplementary Figure 1
% ====================================

lh_similarity = FC_similarity_map(1:40962);
rh_similarity = FC_similarity_map(40963:end);
lh_similarity(isnan(lh_similarity)) = 0;
rh_similarity(isnan(rh_similarity)) = 0;

CBIG_DrawSurfaceMaps(lh_similarity, rh_similarity, 'fsaverage6', 'inflated', 0, 1);
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% ====================================
%% Cleanup
% ====================================

rmpath(fullfile(CBIG_CODE_DIR, 'utilities/matlab/figure_utilities'));
