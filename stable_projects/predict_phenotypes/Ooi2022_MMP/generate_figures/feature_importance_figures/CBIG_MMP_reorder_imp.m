function imp_reordered = CBIG_MMP_reorder_imp(importance, networks)

% CBIG_MMP_reorder_imp(importance, networks)
% Changes the order of the importance vector according to the desired Schaefer400 parcellation order.
% Change the ordering from either the 7 or 17 Yeo networks to the Kong2022 network order.
%
% Input:
% - importance
% A #feature-length vector with importance values. Must be the 79800 vector
% from the Schaefer 400 parcellation.
%
% - networks
% Ordering of old networks. Choose 7 for Yeo7 or 17 for Yeo17 ordering.
% Networks not reordered for any other number.
%
% Output: 
% - imp_reordered
% Reordered importance values.
% 
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% add path to parcellation files
parc_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations');
addpath(fullfile(parc_dir, 'Code'))
label_dir = fullfile(parc_dir, 'FreeSurfer5.3', 'fsaverage6', 'label');
lh_new = fullfile(label_dir, 'lh.Schaefer2018_400Parcels_Kong2022_17Networks_order.annot');
rh_new = fullfile(label_dir, 'rh.Schaefer2018_400Parcels_Kong2022_17Networks_order.annot');

% start reordering vector
imp_mat = CBIG_MMP_FC_vector_2_mat(importance); % reshape into 400x400
% reorder if needed
if networks == 17
    lh_old = fullfile(label_dir, 'lh.Schaefer2018_400Parcels_17Networks_order.annot');
    rh_old = fullfile(label_dir, 'rh.Schaefer2018_400Parcels_17Networks_order.annot');
    idx = CBIG_gwMRF_index_trans_btwn2versions(lh_old, rh_old, lh_new, rh_new);
elseif networks == 7
    lh_old = fullfile(label_dir, 'lh.Schaefer2018_400Parcels_7Networks_order.annot');
    rh_old = fullfile(label_dir, 'rh.Schaefer2018_400Parcels_7Networks_order.annot');
    idx = CBIG_gwMRF_index_trans_btwn2versions(lh_old, rh_old, lh_new, rh_new);
else
    fprintf('Initial network ordering is not Schaefer_Yeo7 or Schaefer_Yeo17, vector not reordered.\n')
end
% extract lower triangle again
imp_mat = imp_mat(idx,idx);
imp_idx = logical(tril(ones(size(imp_mat)), -1));
imp_reordered = imp_mat(imp_idx);

% remove path to parcellation files
rmpath(fullfile(parc_dir, 'Code'))
end