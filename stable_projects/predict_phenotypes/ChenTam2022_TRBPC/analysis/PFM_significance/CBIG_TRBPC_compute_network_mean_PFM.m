function PFM_network = CBIG_TRBPC_compute_network_mean_PFM(PFM_dir, N_score, N_state, outdir)

% PFM_network = CBIG_TRBPC_compute_network_mean_PFM(PFM_dir, N_score, N_state, outdir)
%
% This function computes the predictive feature-matrices averaged within
% each network blocks. It only works for functional connectivity defined on
% 419 (400 cortical and 19 subcortical) Schaefer parcelations
%
% Inputs:
%   - PFM_dir
%     Directory where the predictive-feature matrices are saved
%
%   - N_score
%     Number of behavior scores
%
%   - N_state
%     Number of brain states
%
%   - outdir (optional)
%     Output directory
%
% Outputs:
%   - PFM_network
%     Predictive-feature matrices averaged within each network block
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

N_roi = 419;
N_edge = N_roi*(N_roi-1)/2;

PFM_mean = zeros(N_edge*N_state,N_score);
for i = 1:N_score
    load(fullfile(PFM_dir,['/PFM_score' num2str(i) '_all_folds.mat']),'PFM_all_folds');
    PFM_mean(:,i) = mean(PFM_all_folds,2);
end

N_network = 18;
N_network_pair = N_network*(N_network+1)/2;
PFM_network = zeros(N_network_pair*N_state,N_score);
for i = 1:N_state
    for j = 1:N_score
        PFM_curr_state = PFM_mean((i-1)*N_edge+1:i*N_edge,j);
        PFM_network((i-1)*N_network_pair+1:i*N_network_pair,j) = ...
            CBIG_TRBPC_compute_FC_average_within_network(CBIG_TRBPC_FC_vector2mat(PFM_curr_state));
    end
end

if exist('outdir','var')
    save(fullfile(outdir, 'PFM_network_mean.mat'), 'PFM_network');
end