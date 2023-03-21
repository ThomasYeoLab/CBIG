function [index, parcel_names] = CBIG_ReorderParcelIndex(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot)

% [index, parcel_names] = CBIG_ReorderParcelIndex(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot)
%
% We now have different versions of Schaefer2018 and Yan2023 parcellations.
% For the same resolution, the parcel boundaries are the same but the  
% parcels might be assgined to different networks so the parcel labels
% might be different. Some derived results only depend on the parcel
% boundaries like FC matrices so no need to recompute when switching
% between different versions. 
% This function generates an indexing vector to change the order of FC 
% matrix when switching between different orderings of the same 
% parcellation.
% 
% You can get your new FC matrix by using the output index: 
%       FC_new = FC_old(index, index);
%
% Note that this assumes your FC matrix has the same size as Schaefer
% or Yan parcellation. For example, for 400 Parcels, your FC matrix should
% be 400x400. If you have more ROIs, please manually change the index to 
% match your matrix. 
% For example, you are using Schaefer400 with additional 19 subcortical 
% ROIs, please use:
%       index = [index; (401:419)'];
%       FC_new = FC_old(index, index);
% to keep 19 subcortical ROIs when changing to different versions. 
%
% Input:
%      -lh_old_annot: 
%       Full path of the left hemisphere annot file, the version you used 
%       to compute FC matrix.
%
%      -rh_old_annot:
%       Full path of the right hemisphere annot file, the version you used 
%       to compute FC matrix.
% 
%      -lh_new_annot:
%       Full path of the left hemisphere annot file, the version you want
%       to change to.
%
%      -rh_new_annot:
%       Full path of the right hemisphere annot file, the version you want
%       to change to.
%
% Output:
%      -index:
%       P x 1 vector, Index for following changes. P is the total number of
%       the parcels.
%
%      -parcel_names:
%       Table with new_parcel_name and old_parcel_name, showing the
%       corresponding parcel names between the old version and new version.
% 
% Example:
% [index, parcel_names] = CBIG_ReorderParcelIndex(lh_old_annot, rh_old_annot, lh_new_annot, rh_new_annot)
% FC_new = FC_old(index, index);
% 
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Read annotation files
[~, lh_label_old, lh_table_old] = read_annotation(lh_old_annot);
[~, rh_label_old, rh_table_old] = read_annotation(rh_old_annot);
[~, lh_label_new, lh_table_new] = read_annotation(lh_new_annot);
[~, rh_label_new, rh_table_new] = read_annotation(rh_new_annot);

% Check whether the parcel numbers are the same
if(length(unique(lh_label_old)) ~= length(unique(lh_label_new)))
    error('LH annot files have different parcel numbers.');
end
if(length(unique(rh_label_old)) ~= length(unique(rh_label_new)))
    error('RH annot files have different parcel numbers.')
end

% Extract parcel values
lh_old_label_list = lh_table_old.table(2:end, 5);
rh_old_label_list = rh_table_old.table(2:end, 5);
lh_new_label_list = lh_table_new.table(2:end, 5);
rh_new_label_list = rh_table_new.table(2:end, 5);

lh_index = zeros(size(lh_new_label_list));
rh_index = zeros(size(rh_new_label_list));

for i = 1:length(lh_new_label_list)
    % Find the old label value of the corresponding parcel.
    old_label = unique(lh_label_old(lh_label_new == lh_new_label_list(i)));
    if(length(old_label) ~= 1)
        error('LH annot files have different parcel boundaries.');
    end
    % Find the index of the old parcel
    lh_index(i) = find(lh_old_label_list == old_label);
end
for i = 1:length(rh_new_label_list)
    % Find the old label value of the corresponding parcel.
    old_label = unique(rh_label_old(rh_label_new == rh_new_label_list(i)));
    if(length(old_label) ~= 1)
        error('RH annot files have different parcel boundaries.');
    end
    % Find the index of the old parcel
    rh_index(i) = find(rh_old_label_list == old_label);
end

% Combine left and right
index = [lh_index; rh_index + length(lh_new_label_list)];

if(isequal(index', 1:length(index)))
    disp('The ordering is not changed.');
end

% Create table for corresponding parcel names
old_parcel_name = [lh_table_old.struct_names(2:end); rh_table_old.struct_names(2:end)];
new_parcel_name = [lh_table_new.struct_names(2:end); rh_table_new.struct_names(2:end)];
old_parcel_name = old_parcel_name(index);

parcel_names = table(new_parcel_name, old_parcel_name);

end
