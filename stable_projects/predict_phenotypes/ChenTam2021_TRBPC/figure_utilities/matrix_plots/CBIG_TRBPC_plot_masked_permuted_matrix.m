 function CBIG_TRBPC_plot_masked_permuted_matrix(path_sig_edges, path_all_edges, path_out, lim_maps, cluster_type)
 % CBIG_TRPBC_plot_masked_permuted_matrix(path_sig_edges, path_all_edges, path_out, lim_maps)
 %
 % This function masks out the significant edges (after permutation testing)
 %
 % Required inputs:
 % - path_sig_edges : a filepath to a .mat file containing two 4D matrices of
 %      size (419 x 419 x N conditions x M behaviors), called 
 %      `hyp_driven_mask` and `data_driven_mask`. Each 419 x 419 matrix
 %      contains a binary mask of which edges are significant
 % - path_all_edges : a filepath to a .mat file containing a structure
 %      `mean_stack_clus` that contains the average map for each fmri
 %      condition and behavioral domain
 % - path_out : a filepath to store the outputs
 % - lim_maps : a vector of the lower and upper limits for the output
 %      matrices. If set to [], the default values will be [-2 2]
 % - cluster_type : a string, either 'hypothesis' or 'datadriven', to apply
 %      the masks for either the hypothesis-driven or data-driven behavioral
 %      clusters
 %
 % Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

 
%% set the variables
load(path_sig_edges)
load(path_all_edges)

if strcmp(cluster_type, 'hypothesis')
    sig_edge_task_cat = hyp_driven_mask;
elseif strcmp(cluster_type, 'datadriven')
    sig_edge_task_cat = data_driven_mask;
end

if isempty(lim_maps)
    lim_maps = [-2 2];
end

mkdir(path_out)

fmri_conds = {'rs','mid','nback','sst'};
behav_cat = {'cog','ment_health','pers'};

mean_matrix_cat_sig = [];

% for each behavioral category
for bb = 1:length(behav_cat)
    behav = behav_cat{bb};
    % for each fmri condition
    for ff = 1:length(fmri_conds)
        fmri = fmri_conds{ff};
        % get the matrix
        mat_tmp = mean_stack_clus.(behav).(fmri);
        % get the lower triangle
        tmp_mat_tril = tril(mat_tmp,-1);
        % vectorize
        tmp_vec = nonzeros(tmp_mat_tril(:));
        % divide by the std for visualization
        vec_std = tmp_vec/std(tmp_vec);
        % transform back to matrix
        mat_std = CBIG_TRBPC_FC_vector_2_mat(vec_std);
        % mask out non-significant edges
        mask = sig_edge_task_cat(:,:,ff,bb);
        mat_masked = mat_std.*mask;
        mat_masked_raw = mat_tmp.*mask;
        % plot the one divided by std for visualization
        CBIG_TRBPC_plot_fc_matrix(mat_masked,lim_maps,'hot_cold');
        fname = strcat(path_out, filesep, fmri, '_', behav, '_masked_sig');
        saveas(gcf, fname, 'svg')
        % save the matrix (not divided by std)
        mean_matrix_cat_sig_raw.(behav).(fmri) = mat_masked_raw;
        % save the matrix (divided by std)
        mean_matrix_cat_sig_std.(behav).(fmri) = mat_masked;
    end
end

% save the masked matrices
mat_out = [path_out filesep '/matrices_cat_masked_sig.mat'];
save(mat_out,'mean_matrix_cat_sig_std', 'mean_matrix_cat_sig_raw');



