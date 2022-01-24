function CBIG_TRBPC_plot_avg_relevance(path_fc,path_list_var,...
    list_cog,list_pers,list_ment_health,scalelim,path_out)
% CBIG_TRBPC_plot_avg_relevance(path_fc,path_list_var,...
%     list_cog,list_pers,list_ment_health,scalelim,path_out)
% 
% This function will average together FC relevance values and plot them to produce
% the predictive feature matrices
% 
% Required inputs:
% - path_fc: a string for a path to a .mat file that contains a 
%         structure `struct_fc_vec` with the following fields:
%         `rs`, `mid`, `sst`, `nback`, where each field contains a matrix of
%         size 87571 x b, where b is the number of predicted behaviors
% - path_list_var: a string for a path to a .txt file that contains a list
%         of the predicted behaviors
% - list_cog: a cell array of strings for the cognitive cluster, 
%         where each string is one predicted behavior from `path_list_var`
% - list_pers: a cell array of strings for the personality cluster,
%         where each string is one predicted behavior from `path_list_var`
% - list_ment_health: a cell array of strings for the mental health cluster,
%         where each string is one predicted behavior from `path_list_var`
% - scalelim: a vector containing two numbers specifying the minimum and maximum
%         values for the plots. If left empty, the default is [-2 2]
% - path_out: a string for a path to the output director
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% add required paths to functions
addpath(fullfile(getenv('CBIG_CODE_DIR'),...
    'stable_projects','disorder_subtypes','Tang2020_ASDFactors',...
    'step3_analyses','utilities'));

if ~exist(path_out), mkdir(path_out); end

load(path_fc)

if isempty(scalelim)
    scalelim = [-2 2];
end

% get the names of the behaviours from path_list_var
% read text file line by line and write it into cell_array
ll = 0;
fid = fopen(path_list_var);
while (~feof(fid))
    ll = ll + 1;
    list_behav_orig{ll} = fgetl(fid);
end
fclose(fid);

% matrix of size n by m to store the average FC vectors across 
% fMRI condition (n=4) and cluster (m=3)
vec_stacked = zeros(87571,12);
vv = 0;

% structure to store results of average matrices
mean_stack_clus = [];

% save the new lists into a structure
mean_stack_clus.lists.cog = list_cog;
mean_stack_clus.lists.pers = list_pers;
mean_stack_clus.lists.ment_health = list_ment_health;

% the fmri conditions
conds = {'rs','mid','nback','sst'};
% behavioral clusters
behavs = {'cog','pers','ment_health'};

for bb = 1:length(behavs)
    behav = behavs{bb};
    %% use the clustered list to get the indices from the original full list
    % of bheaviors
    tmp_list = mean_stack_clus.lists.(behav);
    ind_behav = [];
    for ii = 1:length(tmp_list)
        idx = find(strcmp(list_behav_orig,tmp_list{ii}));
        ind_behav(ii) = idx;
    end
    %% for each fmri condition, grab the FC values from tmp_list and plot 
    % the average
    for cc = 1:length(conds)
        cond = conds{cc};
        % get all the vectors for the behaviors in list_cog
        stack = struct_fc_vec.(cond)(:,ind_behav);
        % average them together
        mean_stack = mean(stack,2);
        vv = vv + 1;
        % store the vector of the average
        vec_stacked(:,vv) = mean_stack;
        % divide by the std
        vec_std = mean_stack/std(mean_stack);
        % transform vector to matrix
        mean_stack_mat = CBIG_TRBPC_FC_vector_2_mat(mean_stack);
        mat_std = CBIG_TRBPC_FC_vector_2_mat(vec_std);
        % plot matrix divided by std
        CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(mat_std,scalelim)
        fname = strcat(path_out, filesep, cond, '_', behav, '_clus_sig');
        saveas(gcf,fname,'svg')
        % save rearranged matrix (not divided by std)
        mean_stack_mat = CBIG_TRBPC_rearrange_matrix(mean_stack_mat, 17);
        % store matrix
        mean_stack_clus.(behav).(cond) = mean_stack_mat;
    end
end

%% save mean_stack_clus
fname = strcat(path_out, filesep, 'mean_stack_clus.mat');
save(fname, 'mean_stack_clus')

%% save vec_stacked
fname_vec = strcat(path_out, filesep, 'mean_vec_clus_fmri.mat');
save(fname_vec, 'vec_stacked')

%% remove paths
rmpath(fullfile(getenv('CBIG_CODE_DIR'),...
    'stable_projects','disorder_subtypes','Tang2020_ASDFactors',...
    'step3_analyses','utilities'));
