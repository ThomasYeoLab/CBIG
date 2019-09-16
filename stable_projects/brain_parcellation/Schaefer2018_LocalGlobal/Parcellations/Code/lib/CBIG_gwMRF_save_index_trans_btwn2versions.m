function index_changes = CBIG_gwMRF_save_index_trans_btwn2versions(new_dir)

% index_changes = CBIG_gwMRF_save_index_trans_btwn2versions(new_dir)
%
% This function generates a cell array to save the index changes for all
% Schaefer2018 parcellations.
%
% Input:
%      -new_dir: 
%       Path of the new parcellation folder
%
% Output:
%      -index_changes: 
%       20 x 3 Cell array. The first column is a string showing the
%       resolution of Schaefer2018 parcellations. The second column is a
%       P x 1 vector indicating index changes for each resolution. The
%       third column is a string, if the parcellation order is not changed,
%       the corresponding string will be 'No changes', otherwise it will be
%       'Changed'. 
% 
% Example:
% index_changes = CBIG_gwMRF_save_index_trans_btwn2versions(new_dir)
% 
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'Code'));
output_dir = fullfile(new_dir, 'Updates');
if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end

index_changes = cell(20, 3);
i = 0;
for k = 7:10:17
    for p = 100:100:1000
        i = i + 1;
        lh_old_annot = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
            'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        rh_old_annot = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
            'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        lh_new_annot = fullfile(new_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['lh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        rh_new_annot = fullfile(new_dir, 'FreeSurfer5.3', 'fsaverage', 'label', ...
            ['rh.Schaefer2018_' num2str(p) 'Parcels_' num2str(k) 'Networks_order.annot']);
        index_changes{i, 1} = {[num2str(p) 'Parcels_' num2str(k) 'Networks']};
        index_changes{i, 2} = CBIG_gwMRF_index_trans_btwn2versions(lh_old_annot, rh_old_annot, ...
            lh_new_annot, rh_new_annot);
        if(isequal(index_changes{i, 2}', 1:p))
            index_changes{i, 3} = 'No changes';
        else
            index_changes{i, 3} = 'Changed';
        end
    end
end
save(fullfile(output_dir, 'index_updates.mat'), 'index_changes');

rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'Code'));
end
