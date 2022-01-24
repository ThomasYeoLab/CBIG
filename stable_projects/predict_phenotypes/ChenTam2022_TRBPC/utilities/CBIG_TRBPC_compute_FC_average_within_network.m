function FC_net = CBIG_TRBPC_compute_FC_average_within_network(FC)

% FC_net = CBIG_TRBPC_compute_FC_average_within_network(FC)
%
% This fucntion computes the functional connectivity average within 18*18
% network blocks
%
% Inputs:
%   - FC
%     A 419*419 FC matrix. The FC must be an output of CBIG_fMRI_Preproc2016
%     pipeline
%
% Outputs:
%   -FC_mean_net
%    FC averaged within each 18*18 network pair blocks. We only take the
%    lower triangular and diagonal entries of this 18*18 matrix so it is a vector of lenght 18*19/2
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

N_network = 18;
N_roi = 419;
load(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', 'ChenTam2022_TRBPC', ...
    'utilities', 'reorder_FC_mat.mat'));
FC_reorder = FC(index_to_reorder,index_to_reorder);

network_18 = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'Limbic_TempPole';...
    'Limbic_OFC'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; ...
    'VisPeri'; 'VisCent'; 'Subcortical'};

index_network = cell(N_network,1);
nodes = 1:N_roi;
for i = 1:N_network-1
    curr_net = network_18{i};
    index_network{i} = nodes(contains(label_after_reorder,curr_net));
end
index_network{18} = 401:419;

FC_mean = zeros(N_network,N_network);
for i = 1:N_network
    for j = 1:N_network
        FC_mean(i,j) = ...
            mean(mean(FC_reorder(index_network{i},index_network{j})));
    end
end

ind = tril(ones(18,18))==1;
FC_net = FC_mean(ind);
