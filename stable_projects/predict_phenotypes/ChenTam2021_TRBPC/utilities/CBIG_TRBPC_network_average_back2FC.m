function FC = CBIG_TRBPC_network_average_back2FC(FC_network)

% FC = CBIG_TRBPC_network_average_back2FC(FC_network)
%
% This fucntion maps 18*18 functional connectivity (FC) network averages back to 419*419 FC matrix
%
% Inputs:
%   - FC_network
%     A 18*18 matrix. Each entry in this matrix represents the mean FC
%     value of a netwok block
%
% Outputs:
%   - FC
%     A 419*419 FC matrix. FC value of edges in each network block in this matrix equals the
%     mean FC value of that block in FC_network
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

N_network = 18;
N_roi = 419;
load reorder_FC_mat.mat

network_18 = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; 'Limbic_TempPole';...
    'Limbic_OFC'; 'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA';...
    'VisPeri'; 'VisCent'; 'Subcortical'};

index_network = cell(N_network,1);
nodes = 1:N_roi;
for i = 1:N_network-1
    curr_net = network_18{i};
    index_network{i} = nodes(contains(label_after_reorder,curr_net));
end
index_network{18} = 401:419;

FC = zeros(N_roi,N_roi);
for i = 1:N_network
    for j = 1:N_network
        FC(index_network{i},index_network{j}) = FC_network(i,j);
    end
end