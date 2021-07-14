function CBIG_TRBPC_aggregate_relevance(path_data,path_out)
% CBIG_TRBPC_aggregate_relevance(path_data,path_out)
%
% this function will aggregate the weights across the folds for each 
% predicted behavior and fMRI condition
%
% required inputs:
% - path_data: a string for a path to a directory that contains 36 .mat
%       files (1 for each of the 36 behaviors) with the naming convention 
%       'activation_score_X_120folds.mat' where X is a number representing 
%       a given behavior
% - path_out: a string for a path to a directory for the outputs
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if ~exist(path_out), mkdir(path_out); end

conds = {'rs','mid','nback','sst'};
conds_ind = {1:87571, 87572:175142, 175143:262713, 262714:350284};

struct_fc_vec = [];

for cc = 1:length(conds)
    cond = conds{cc};
    struct_fc_vec.(cond) = zeros(87571,36);
end

stack = zeros(350284,36);

% for each behavioural score
for bb = 1:36
    path_mat = [path_data filesep 'PFM_score' num2str(bb) '_all_folds.mat'];
    load(path_mat)
    struct_mat = [];
    % average across folds/seeds to get vector
    stack_fc_vec = mean(PFM_all_folds,2);
    stack(:,bb) = stack_fc_vec;
    % for each fmri condition
    for cc = 1:length(conds)
        cond = conds{cc};
        cond_ind = conds_ind{cc};
        score_name = strcat('score',num2str(bb));
        % take one condition
        struct_mat.(cond) = PFM_all_folds(cond_ind,:);
        % average across folds/seeds to get vector
        vec_fc = mean(struct_mat.(cond),2);
        % save the vector
        fprintf(strcat("saving ", cond, ' score', num2str(bb), "\n"))
        struct_fc_vec.(cond)(:,bb) = vec_fc;
    end
end

% save vectors
mat_out = [path_out filesep 'relevance_vectors.mat'];
save(mat_out,'struct_fc_vec')

mat_out = [path_out filesep 'stacked_relevance_vectors.mat'];
save(mat_out,'stack')


