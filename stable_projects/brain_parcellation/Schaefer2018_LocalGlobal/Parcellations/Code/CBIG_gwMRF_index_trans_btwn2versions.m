function index = CBIG_gwMRF_index_trans_btwn2versions(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot)

% index = CBIG_gwMRF_index_trans_btwn2versions(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot)
%
% We now have different versions of Schaefer2018 parcellations. If you want
% to update from v0.8.1-Schaefer2018_LocalGlobal and earlier versions to
% v0.14.3-Update_Yeo2011_Schaefer2018_labelname and future versions, you 
% can use this function to revise your previous calculations, for example, 
% FC matrix.
% 
% This function generates an indexing vector to change the order of FC 
% matrix when you switch from the old annot files to the new annot files. 
% 
% You can get your new FC matrix by using the output index: 
%       FC_new = FC_old(index, index);
%
% Note that this assumes your FC matrix has the same size as Schaefer
% parcellation. For example, for 400 Parcels, your FC matrix should be
% 400x400. If you have more ROIs, please manually change the index to match
% your matrix. 
%
% Input:
%      -lh_old_annot: 
%       Full path of the left hemisphere annot file, the version you use to
%       compute FC matrix.
%
%      -rh_old_annot:
%       Full path of the right hemisphere annot file, the version you use 
%       to compute FC matrix.
% 
%      -lh_new_annot:
%       Full path of the left hemisphere annot file, the version you want
%       to update.
%
%      -rh_new_annot:
%       Full path of the right hemisphere annot file, the version you want
%       to update.
%
% Output:
%      -index:
%       P x 1 vector, Index for following changes. P is the total number of
%       the parcels . 
% 
% Example:
% index = CBIG_gwMRF_index_trans_btwn2versions(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot)
% 
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[~, lh_label_old, lh_table_old] = read_annotation(lh_old_annot);
[~, rh_label_old, rh_table_old] = read_annotation(rh_old_annot);
[~, lh_label_new, lh_table_new] = read_annotation(lh_new_annot);
[~, rh_label_new, rh_table_new] = read_annotation(rh_new_annot);
lh_list = lh_table_new.table(2:end, 5);
rh_list = rh_table_new.table(2:end, 5);
lh_index = zeros(size(lh_list));
rh_index = zeros(size(rh_list));

for i = 1:length(lh_list)
    lh_index(i) = find(lh_table_old.table(:, 5) == mode(lh_label_old(lh_label_new == lh_list(i))));
end
for i = 1:length(rh_list)
    rh_index(i) = find(rh_table_old.table(:, 5) == mode(rh_label_old(rh_label_new == rh_list(i))));
end

index = [lh_index; rh_index + length(lh_list)] - 1;

if(isequal(index', 1:length(index)))
    disp('The ordering is not changed.');
else

end
