function CBIG_TRBPC_supp_relevance_ind(path_fc, lims, opt_std, path_out)
% CBIG_TRBPC_supp_relevance_ind(path_fc, lims, opt_std, path_out)
%
% This function will put all the relevance maps for the scores together in one plot per
% fmri condition
%
% Required inputs:
% - path_fc: a string for a path to a .mat file that contains a structure
%         called `struct_fc_vec` which contains the following subfields
%         `rs`, `mid`, `sst`,`nback`, which each have vectors of the relevance
%         values for each edge of each predicted behavior
% - lims: a vector of two numerical values specifying the minimum and maximum
%         values in the plots, e.g. [-1 1]
% - opt_std: a binary variable, 0 or 1. If set to 1, the plots will be divided by
%         the standard deviation. If set to 0, the raw values will be shown.
% - path_out: a string to specify the output directory
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% add required paths
addpath(fullfile(getenv('CBIG_CODE_DIR'),'external_packages','matlab',...
    'non_default_packages','figure_utilities'));

% output directory
if ~exist(path_out), mkdir(path_out); end

% load data
load(path_fc)
% Load colormap
load corr_mat_colorscale.mat;

% list of behaviors in original order
list_behav_orig = {'anxious depressed','withdrawn depressed','somatic complaints',...
    'social problems','thought problems','attention problems','rulebreaking behavior',...
    'aggressive behavior','vocabulary','attention','working memory',...
    'executive function','processing speed','episodic memory','reading',...
    'fluid cognition','crystallized cognition','overall cognition',...
    'negative urgency','lack of planning','sensation seeking','positive urgency',...
    'lack perseverance','behavioral inhibition','reward responsiveness','drive',...
    'fun seeking','total psychosis symptoms','psychosis severity',...
    'mania severity','short delay recall','long delay recall',...
    'fluid intelligence','visuospatial accuracy','visuospatial reaction time',...
    'visuospatial efficiency'};

% re-sorted list of behaviors (cognition, personality, mental health scales)
list_behav_new = {'vocabulary','attention','working memory','executive function',...
    'processing speed','episodic memory','reading','fluid cognition',...
    'crystallized cognition','overall cognition','short delay recall',...
    'long delay recall','fluid intelligence','visuospatial accuracy',...
    'visuospatial reaction time','visuospatial efficiency','negative urgency',...
    'lack of planning','sensation seeking','positive urgency','lack perseverance',...
    'behavioral inhibition','reward responsiveness','drive','fun seeking',...
    'anxious depressed','somatic complaints','social problems','thought problems',...
    'attention problems','total psychosis symptoms','psychosis severity'};

%% group behaviors to prep for a color scheme
list_nih = {'vocabulary','attention','working memory','executive function',...
    'processing speed','episodic memory','reading','fluid cognition',...
    'crystallized cognition','overall cognition'};
list_recall = {'short delay recall','long delay recall'};
list_fi = {'fluid intelligence'};
list_vis = {'visuospatial accuracy','visuospatial reaction time',...
    'visuospatial efficiency'};
list_upps = {'negative urgency','lack of planning',...
    'sensation seeking','positive urgency','lack perseverance'};
list_bisbas = {'behavioral inhibition','reward responsiveness','drive','fun seeking'};
list_cbcl = {'anxious depressed','somatic complaints','social problems',...
    'thought problems','attention problems'};
list_pps = {'total psychosis symptoms','psychosis severity'};

%% make the plots
fmris = {'rs','mid','nback','sst'};

for ff = 1:length(fmris)
    close all
    fmri = fmris{ff};
    %figure
    for bb = 1:length(list_behav_new)
        behav = list_behav_new{bb};
        % get index from old list
        idx = find(strcmp(list_behav_orig,behav));
        % get the vector
        tmp_vec = struct_fc_vec.(fmri)(:,idx);
        % if user wants to divide by the std, do so
        if opt_std == 1
            tmp_vec = tmp_vec/std(tmp_vec);
        end
        % transform to the matrix
        tmp_mat = CBIG_TRBPC_FC_vector_2_mat(tmp_vec);
        % rearrange
        tmp_mat = CBIG_TRBPC_rearrange_matrix(tmp_mat,17);
        % plot
        subplot_tight(7, 5, bb, [0.02, 0.02])
        imshow(tmp_mat,lims)
        colormap(gca,rbmap2)
        if any(find(strcmp(list_nih,behav)))
            col_vec = [1     0     0];
        elseif any(find(strcmp(list_recall,behav)))
            col_vec = [0.6980    0.1333    0.1333];
        elseif any(find(strcmp(list_fi,behav)))
            col_vec = [1.0000    0.2706         0];
        elseif any(find(strcmp(list_vis,behav)))
            col_vec = [0.9412    0.5020    0.5020];
        elseif any(find(strcmp(list_upps,behav)))
            col_vec = [0.1922    0.1922    0.1922];
        elseif any(find(strcmp(list_bisbas,behav)))
            col_vec = [0.5882    0.5451    0.5451];
        elseif any(find(strcmp(list_cbcl,behav)))
            col_vec = [0     0     1];
        elseif any(find(strcmp(list_pps,behav)))
            col_vec = [0.3922    0.5843    0.9294];
        end
        title(behav, 'FontWeight', 'Normal', 'FontName', 'Arial', 'FontSize', 8, 'Color', col_vec)
    end
    fname = [path_out filesep 'sig_all_' fmri];
    print('-fillpage',fname,'-dpdf')
    set(gcf, 'PaperUnits', 'inches')
    set(gcf, 'PaperPosition',[0 0 11 8.5])
    saveas(gcf,fname,'svg');
    close all
end

%% remove path
rmpath(fullfile(getenv('CBIG_CODE_DIR'),'external_packages','matlab',...
    'non_default_packages','figure_utilities'));
        
    