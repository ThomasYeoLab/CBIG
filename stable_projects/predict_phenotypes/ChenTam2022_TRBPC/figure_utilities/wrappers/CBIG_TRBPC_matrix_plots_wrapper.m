function CBIG_TRBPC_matrix_plots_wrapper(data_dir, top_outdir, hypothesis_ind, datadriven_ind, significant_edges)
% CBIG_TRBPC_matrix_plots_wrapper(data_dir, top_outdir, hypothesis_ind, datadriven_ind)
% 
% This function will generate all matrix-related figures for the TRBPC project
%
% Inputs:
%   - data_dir
%     Directory where regression and PFM results are stored
%
%   - top_outdir
%     Top level output directory
%
%   - hypothesis_ind
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for hypothesis-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
%   - datadriven_ind
%     A structure of 3 fields:
%         .cog: a vector containing the indices of cognitive measures
%               for data-driven clusters
%         .pers: a vector containing the indices of personality measures
%         .mh: a vector containing the indices of mental health measures
%
%   - significant_edges
%   A .mat file containing the mask for significant edges. The file is an
%   output of `CBIG_TRBPC_compute_PFM_significance_mask_wrapper.m`
%
% Outputs:
%   Output figures and .mat files will be stored in top_outdir
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% add required paths
addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects',...
    'predict_phenotypes','ChenTam2022_TRBPC','figure_utilities')));

in_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects',...
    'predict_phenotypes','ChenTam2022_TRBPC','figure_utilities','input');

%% make output dir
if ~exist(top_outdir); mkdir(top_outdir); end

%% step 1: aggregate the weights across the folds for each predicted 
% behavior and fMRI condition
in = fullfile(data_dir,'PFM','KRR','allFC');
out = top_outdir;

CBIG_TRBPC_aggregate_relevance(in, out);
clear in out

%% step 2: average together FC relevance values and plot them

cluster_types = {'datadriven', 'hypothesis'};

in.path_fc = [top_outdir filesep 'relevance_vectors.mat'];
in.path_list_var = [in_dir filesep 'variables_to_predict.txt'];
list_var = CBIG_text2cell(in.path_list_var);

for cc = 1:length(cluster_types)
    clus = cluster_types{cc};
    if strcmp(clus, 'hypothesis')
        in.list_cog = list_var(hypothesis_ind.cog);
        in.list_pers = list_var(hypothesis_ind.pers);
        in.list_ment_health = list_var(hypothesis_ind.mh);
    elseif strcmp(clus, 'datadriven')
        in.list_cog = list_var(datadriven_ind.cog);
        in.list_pers = list_var(datadriven_ind.pers);
        in.list_ment_health = list_var(datadriven_ind.mh);
    end
    in.scalelim = [-2 2];
    out.path_out = [top_outdir filesep clus];
    
    CBIG_TRBPC_plot_avg_relevance(in.path_fc, in.path_list_var,...
        in.list_cog, in.list_pers, in.list_ment_health, in.scalelim, out.path_out)
    
    close all
end
clear in out

%% step 3: mask out the edges that are not significant

for cc = 1:length(cluster_types)
    clus = cluster_types{cc};

    in.path_sig_edges = significant_edges;
    in.path_all_edges = [top_outdir filesep clus filesep 'mean_stack_clus.mat'];
    out.path_out = [top_outdir filesep clus];
    opt.lim_maps = [-1.5 1.5];
    
    if strcmp(clus, 'hypothesis')
        opt.cluster_type = 'hypothesis';
    elseif strcmp(clus, 'datadriven')
        opt.cluster_type = 'datadriven';
    end
    
    CBIG_TRBPC_plot_masked_permuted_matrix(in.path_sig_edges, in.path_all_edges,...
        out.path_out, opt.lim_maps, opt.cluster_type)
    
    close all
    clear in out opt
end

%% step 4: do a conjunction analysis and plot matrices of the conjunction

for cc = 1:length(cluster_types)
    clus = cluster_types{cc};
    
    in.path_data = [top_outdir filesep clus];
    opt.lim_maps = [-1.5 1.5];
    
    CBIG_TRBPC_edges_conjunction(in.path_data, opt.lim_maps)
    
    close all
    clear in opt
end

%% step 5: put the significant positive and significant negative blocks 
% together into one figure for the relevance matrix and chord diagrams

for cc = 1:length(cluster_types)
    clus = cluster_types{cc};

    in.path_data = [top_outdir filesep clus filesep 'conjunction'];
    opt.lim_maps = [-1.5 1.5];
    
    CBIG_TRBPC_chord_matrix_merge_pos_neg(in.path_data, opt.lim_maps)
    
    close all
    clear in opt
end

%% step 6: generate .annot files for cortical surface projection

for cc = 1:length(cluster_types)
    clus = cluster_types{cc};
    
    in.path_data = [top_outdir filesep clus filesep 'conjunction' filesep, ...
        'surfaces/vec_relevance_sum_conjunction.mat'];
    out.path_out = [top_outdir filesep clus filesep 'conjunction' filesep, ...
        'surfaces'];
    
    CBIG_TRBPC_cortical_surface_projection(in.path_data, out.path_out)
    clear in out
end

%% step 7: generate subcortical heatmaps

for cc = 1:length(cluster_types)
    clus = cluster_types{cc};

    in.path_data = [top_outdir filesep clus filesep 'conjunction/surfaces'];
    
    % set the value limits for the figures
    load([in.path_data filesep 'vec_relevance_sum_conjunction.mat'])
    behs = {'cog','pers','ment_health'};
    for bb = 1:length(behs)
        beh = behs{bb};
        % grab the vectors
        tmp.vec_pos = vec_relevance_conjunction.pos.(beh);
        tmp.vec_neg = vec_relevance_conjunction.neg.(beh);
        % get the 95th percentile of both
        tmp.max_pos = prctile(tmp.vec_pos,95);
        tmp.max_neg = prctile(tmp.vec_neg,95);
        % make the limit the smaller 95th percentile
        tmp.lim_p_max = min([tmp.max_pos tmp.max_neg]);
        tmp.tmp_str = strcat('lim_', beh);
        opt.(tmp.tmp_str) = [0 tmp.lim_p_max];
    end
    
    CBIG_TRBPC_subcortical_heatmap(in.path_data,...
        opt.lim_cog, opt.lim_pers, opt.lim_ment_health)
    close all
    clear vec_relevance_conjunction
    clear in opt
end

%% step 8: supplemental PFMs for each behavior and fMRI state

in.path_fc = [top_outdir filesep filesep 'relevance_vectors.mat'];
opt.lims = [-2 2];
opt.opt_std = 1;
out.path_out = top_outdir;

CBIG_TRBPC_supp_relevance_ind(in.path_fc, opt.lims, opt.opt_std, out.path_out)
close all

%% rm required paths
rmpath(genpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects',...
    'predict_phenotypes','ChenTam2022_TRBPC','figure_utilities')));

