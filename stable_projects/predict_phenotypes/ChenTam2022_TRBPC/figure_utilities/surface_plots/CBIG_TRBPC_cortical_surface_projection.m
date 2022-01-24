function CBIG_TRBPC_cortical_surface_projection(path_data,path_out)
% CBIG_TRBPC_cortical_surface_projection(path_data,path_out)
%
% This function will create files to generate cortical surface figures
%
% Required inputs:
% - path_data : a path to a .mat file that contains a structure named
%              `vec_relevance_conjunction`
% - path_out : a path to save results
%
% Written by Angela Tam & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% function starts here

% add required paths
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects',...
    'predict_phenotypes','ChenTam2022_TRBPC','figure_utilities'));

% create output directory
if ~exist(path_out), mkdir(path_out); end

% load the data
load(path_data)

vals = {'pos','neg'};
clus = {'cog','pers','ment_health'};

for cc = 1:length(clus)
    clu = clus{cc};
    % grab the vectors
    vec_pos = vec_relevance_conjunction.pos.(clu);
    vec_neg = vec_relevance_conjunction.neg.(clu);
    % get the 95th percentile of both
    max_pos = prctile(vec_pos,95);
    max_neg = prctile(vec_neg,95);
    % make the limit the smaller 95th percentile
    lim_p_max = min([max_pos max_neg]);
    lim_p = [0 lim_p_max];
    % generate the files
    for vv = 1:length(vals)
        val = vals{vv};
        % grab the vector
        vec =  vec_relevance_conjunction.(val).(clu);
        % put onto surface
        out_dir = strcat(path_out, filesep, val, '_', clu, '/annot');
        % set limits to be 5th to 95th percentiles
        if strcmp(val, 'neg')
            CBIG_TRBPC_Ruby_create_annot_r_overlay_Schaefer400(vec, out_dir, lim_p, 'blue_cyan')
        else
            CBIG_TRBPC_Ruby_create_annot_r_overlay_Schaefer400(vec, out_dir, lim_p, 'red_yell')
        end
    end
end

% remove required paths
rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects',...
    'predict_phenotypes','ChenTam2022_TRBPC','figure_utilities'));