function CBIG_IndCBM_extract_MSHBM_result(project_dir)

% CBIG_IndCBM_extract_MSHBM_result(project_dir)
% This function extract individual cerebral cortical parcellations after
% estimating group priors by step 2 of MSHBM algorithm.
%
% Input:
%     - project_dir: 
%           Params_Final.mat from MSHBM step 2 should be saved under
%           project_dir/priors
%           Extracted individual cerebral cortical parcellations will be
%           saved under project_dir/ind_parcellation
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Load mat file
final = load(fullfile(project_dir, 'priors', 'Params_Final.mat'));
vertex_num = size(final.Params.s_lambda, 1) / 2;
num_clusters = size(final.Params.s_lambda, 2);
sub_num = size(final.Params.s_lambda, 3);

output_dir = fullfile(project_dir, 'ind_parcellation');
if(~exist(output_dir))
    mkdir(output_dir);
end

% Extract individual parcellation
for i = 1:sub_num
    if(isempty(final.Params.s_lambda))
        error('s_lambda is empty. Please use save_all flag when estimating group prior.')
    end
    sub = final.Params.s_lambda(:, :, i);
    [~, labels] = max(sub');
    labels(sum(sub, 2) == 0) = 0; % Label medial wall as 0
    lh_labels = labels(1:vertex_num)';
    rh_labels = labels((vertex_num+1):end)';
    output_file = fullfile(output_dir, ['Ind_parcellation_MSHBM_sub' num2str(i) '.mat']);
    group_file = fullfile(project_dir, 'group', 'group.mat');
    if(exist(group_file, 'file'))
        group = load(group_file);
        if(isfield(group, 'colors')) % Use group color table for individual parcellations
            colors = group.colors;
            save(output_file, 'lh_labels', 'rh_labels', 'colors', 'num_clusters');
        else
            save(output_file, 'lh_labels', 'rh_labels', 'num_clusters');
        end
    else
        save(output_file, 'lh_labels', 'rh_labels', 'num_clusters');
    end
end

end
