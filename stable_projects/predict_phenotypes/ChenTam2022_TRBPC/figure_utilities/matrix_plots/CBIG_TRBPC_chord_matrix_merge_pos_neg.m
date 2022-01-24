function CBIG_TRBPC_chord_matrix_merge_pos_neg(path_data, lim_maps)
% CBIG_TRBPC_chord_matrix_merge_pos_neg(path_data, lim_maps)
% 
% This function puts the average values of the mainly positive and mainly 
% negative blocks together into one csv, and
% also plots to a 419 x 419 matrix for each behavioral cluster
%
% required inputs
% - path_data : a filepath to a directory that contains a .mat file called
%       `mat_relevance_conjunction_masked.mat`; this directory will also be
%       used to save the outputs
% - lim_maps : a vector for the lower and upper limits for the output
%       matrices. If set to empty, the defaults will be [-2 2]
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% function starts here
    
% set the variables
path_out = path_data;
path_mat = [path_data filesep 'mat_relevance_conjunction_masked.mat'];
load(path_mat)
% names of behavioral cluster inside path_mat
clus = {'cog','ment_health','pers'};
vals = {'pos','neg'};
% default limits for output matrices
if isempty(lim_maps)
    lim_maps = [-2 2];
end

for cc = 1:length(clus)
    clu = clus{cc};
    matrices_tmp = zeros(419,419,2);
    nn = 1;
    for vv = 1:length(vals)
        val = vals{vv};
        % get the matrix
        matrices_tmp(:,:,nn) = masked_avg_relevance_conjunction.avg.(val).(clu);
        nn = nn + 1;
    end
    % sum the matrices element wise
    matrix_sum = sum(matrices_tmp,3);
    % plot the matrix
    CBIG_TRBPC_plot_fc_matrix(matrix_sum,lim_maps,'hot_cold')
    fname = strcat(path_out, filesep, 'relevance_masked_conjunction_bothdir_', clu);
    saveas(gcf,fname,'svg')
    % make logical to show the conjunction and plot
    CBIG_TRBPC_plot_fc_matrix(logical(matrix_sum),[0 1],'hot')
    fname = strcat(path_out, filesep, 'conjunction_bothdir_', clu);
    saveas(gcf,fname,'svg')
    % average the values within network blocks for chord diagrams
    avg_rel_blk = CBIG_TRBPC_avg_network_pairs_relevance(matrix_sum);
    % write to csv
    csv_out = [path_out filesep 'chord/avg_bothdir_' clu '.csv'];
    csvwrite(csv_out, avg_rel_blk)
end
    
        
    

