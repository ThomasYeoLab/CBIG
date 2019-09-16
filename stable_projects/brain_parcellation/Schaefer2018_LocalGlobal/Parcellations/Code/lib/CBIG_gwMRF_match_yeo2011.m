function [label_matched, table_matched] = ...
    CBIG_gwMRF_match_yeo2011(alex_annot, yeosplit_annot, hemi, network_name, k, p, reorder_idx)

% [label_matched, table_matched] = ...
%       CBIG_gwMRF_match_yeo2011(alex_annot, yeosplit_annot, hemi, network_name, k, p, reorder_idx)
%
% This function match the Schaefer2018 parcellation with the Yeo2011 split
% components. For each parcel, we first assgin the network, then assgin 
% the component. 
%
% Input:
%      -alex_annot: 
%       Path of Schaefer2018 parcellation annot file.
%
%      -yeosplit_annot:
%       Path of growed Yeo2011 split components annot file. 
% 
%      -hemi:
%       'lh' or 'rh'
%
%      -network_name:
%       Cell arrays where each element contains the name of a network.
%
%      -k: 
%       A Scalar. Total number of networks.
% 
%      -p: 
%       A Scalar. Total number of parcels.
%
%      -reorder_idx:
%       A vector to change the order of the components. 
%       Example: [1:28,30,29,31:57]
%       If not given, then the order of the components will be the same as
%       the order of the color table of Yeo2011 split components annot file. 
%
% Output:
%      -label_matched:
%       Matched labels on fsaverage. 
% 
%      -table_matched:
%       Ordered color table
%
% Example:
% [label_matched, table_matched] = ...
%       CBIG_gwMRF_match_yeo2011(alex_annot, yeosplit_annot, 'lh', network_name, 7, 400)
%       
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[~, label_yeosplit, table_yeosplit] = read_annotation(yeosplit_annot);
[~, label_alex, table_alex] = read_annotation(alex_annot);
label_alex = label_alex';

if(exist('reorder_idx','var'))
    table_yeosplit.table(2:end, :) = table_yeosplit.table(1 + reorder_idx, :);
    table_yeosplit.struct_names(2:end) = table_yeosplit.struct_names(1 + reorder_idx);
end

% Read mesh
mesh7 = CBIG_ReadNCAvgMesh(hemi, 'fsaverage', 'sphere', 'cortex');
mesh5 = CBIG_ReadNCAvgMesh(hemi, 'fsaverage5', 'sphere', 'cortex');
surface = CBIG_ReadNCAvgMesh(hemi, 'fsaverage', 'inflated', 'cortex');

% Upsample to fs7
label_yeosplit = MARS_NNInterpolate_kdTree(mesh7.vertices, mesh5, label_yeosplit');

% Remove medial wall
label_yeosplit(label_yeosplit == 65793) = nan;

label_yeo = label_yeosplit;
for i=2:length(table_yeosplit.struct_names)
    component_index = find(label_yeosplit == table_yeosplit.table(i,5));
    component_name = table_yeosplit.struct_names{i};
    for j = 2:length(network_name)
        if (contains(component_name, network_name{j}))
            network = j;
            break
        end
    end
    label_yeo(component_index) = network;
end

label_max_match = zeros(size(label_alex));
label_alex_list = unique(label_alex);
label_alex_list(label_alex_list == table_yeosplit.table(1,5)) = [];
label_alex_list(label_alex_list == 0) = [];
for i=1:length(label_alex_list)
    network = mode(label_yeo(label_alex == label_alex_list(i)));
    component_labels = label_yeosplit(label_alex == label_alex_list(i));
    component_list = unique(component_labels);
    component_list(isnan(component_list)) = [];
    count = 0;
    ignore_list = zeros(size(component_labels));
    for j = 1:length(component_list)
        index = find(table_yeosplit.table(:, 5) == component_list(j));
        if(~contains(table_yeosplit.struct_names{index}, network_name{network}))
            ignore_list(component_labels == component_list(j)) = 1;
            count = count + 1;
        end
    end
    if(count<length(component_list) && sum(ignore_list)>1)
        component_labels(logical(ignore_list)) = nan;
    end
    if(count == length(component_list))
        warning_index = find(table_alex.table(:,5)==label_alex_list(i));
        warning(['Please check ' hemi ', parcel ' table_alex.struct_names{warning_index} '.']);
    end
    label_max_match(label_alex==label_alex_list(i))=mode(component_labels);
end

% This is a manual change
label_max_match = manual_changes(label_max_match, label_alex, table_yeosplit, table_alex, k, p);

% Assign names to different parcels within each component
[label_matched, table_matched]=...
    alex_reorder_parcel_numbers_across_hemispheres(...
    label_alex, table_yeosplit, label_max_match, surface.vertices(3,:));

end

function label_new = manual_changes(label_matched, label_alex, table_yeosplit, table_alex, k, p)

% label_new = manual_changes(label_matched, label_alex, table_yeosplit, table_alex, k, p)
% 
% This function does some manual changes after assigning components since
% some of the parcels can not be automatically labeled correctly. If the
% user wish to add further manual changes, please refer to the first line
% of the main code and make the necessary amendments.
%
% For example, from version v0.8.1-Schaefer2018_LocalGlobal to version
% v0.14.3-Update_Yeo2011_Schaefer2018_labelname, parcel 
% '7Networks_LH_Default_Temp_3' in the 300-parcel, 7-network parcellation 
% are assigned to '7Networks_LH_SomMot' network. Since we don't want to 
% change its network label, we manually change its assignment into 
% '7Networks_LH_Default_Temp'. 
% 
% Input:
%      -label_matched: 
%       A column vector, automatically assigned labels
% 
%      -label_alex: 
%       A column vector, label of old Schaefer2018 parcellation
% 
%      -table_yeosplit: 
%       Color table of Yeo2011 split components
% 
%      -table_alex: 
%       Color table of old Schaefer2018 parcellation
% 
%      -k: 
%       A Scalar. Total number of networks.
% 
%      -p: 
%       A Scalar. Total number of parcels.
%
% Example:
% label_new = manual_changes(label_matched, label_alex, table_yeosplit, table_alex, 17, 400)

% If you want to add manual changes, please change the following variables:
resolusion = [7, 300]; % Each row corresponds to a single change. 
            % The first column is the network index and the second column is the parcel index
old_parcel{1} = '7Networks_LH_Default_Temp_3'; % old_parcel{i} is the old name of the ith change
new_component{1} = '7Networks_LH_Default_Temp'; % new_component{i} is the new component name of the i-th change

label_new = label_matched;
for i = 1:size(resolusion, 1)
    if(k == resolusion(i, 1) && p == resolusion(i, 2))
        index_alex = 0;
        index_yeo = 0;
        
        % Find the label index of the parcel
        for j = 1:length(table_alex.struct_names)
            if(strcmp(table_alex.struct_names{j}, old_parcel{1}))
                index_alex = j;
                break
            end
        end
        
        % Find the label index of the component
        for j = 1:length(table_yeosplit.struct_names)
            if(strcmp(table_yeosplit.struct_names{j}, new_component{1}))
                index_yeo = j;
                break
            end
        end
        
        % Assign the parcel with the new component label
        if(index_alex && index_yeo)
            parcel = label_alex == table_alex.table(index_alex, 5);
            label_new(parcel) = table_yeosplit.table(index_yeo, 5);
            disp(['Parcel ' table_alex.struct_names{index_alex} ' is manually labeled as ' ...
                table_yeosplit.struct_names{index_yeo} ' for ' num2str(p) ' Parcels. ']);
        end
    end
end

end

function [label_new, table_new] = ...
    alex_reorder_parcel_numbers_across_hemispheres(label_alex, table_yeosplit, label_max_match, coordinates)

% [label_new, table_new] = ...
%     alex_reorder_parcel_numbers_across_hemispheres(label_alex, table_yeosplit, label_max_match, coordinates)
% 
% This function numbers the parcels within the same components according to
% the coordinates and also color the parcels with different colors.
% 
% Input:
%      -label_alex: 
%       A column vector, label of old Schaefer2018 parcellation
% 
%      -table_yeosplit: 
%       Color table of Yeo2011 split components
% 
%      -label_max_match: 
%       Component labels of all the Schaefer Parcellation vertices after
%       matching to Yeo2011 split components labels. Vertices in the same
%       components have the same labels.
% 
%      -coordinates: 
%       Coordinate of each vertex, this is used for numbering.
%
% Example:
% [label_new, table_new] = ...
%     alex_reorder_parcel_numbers_across_hemispheres(label, table, label_max_match, surface.vertices(3,:))

label_new = zeros(size(label_max_match));

% neutral element, this might need to be changed in other settings
table_new.table(1, :) = table_yeosplit.table(1, :); 
table_new.struct_names(1) = table_yeosplit.struct_names(1, :);

current_label = 1;% this will be upcounted
for i = 2:length(table_yeosplit.table(:, 5))
    current_color = table_yeosplit.table(i, 5);
    [current_orig_labels] = unique(label_alex(label_max_match == current_color));
    positions = zeros(1, length(current_orig_labels));
    % Ordering within original network parcel
    for j = 1:length(current_orig_labels)
        positions(j) = mean(coordinates(label_alex == current_orig_labels(j)));
    end
    [sorted_pos, idx] = sort(positions);
    internal_label = 1; % these are (sub) parcels within a network
    off_set_internal = -ceil(length(sorted_pos) / 2); % for collision check
    for j = 1:length(sorted_pos)
        current_label = current_label + 1;
        table_new.struct_names(current_label) = strcat(table_yeosplit.struct_names(i), '_', num2str(internal_label));
        internal_label = internal_label + 1;
        table_new.table(current_label,:) = table_yeosplit.table(i, :);
        
        vec = compute_vec_from_label(internal_label + off_set_internal);
        collision = check_for_collision(vec, table_new, current_label);
        while(collision > 0)
            off_set_internal = off_set_internal + 1;
            vec = compute_vec_from_label(internal_label + off_set_internal);
            collision = check_for_collision(vec, table_new, current_label);
            off_set_internal = check_for_overflow(vec, table_new, current_label, off_set_internal);
        end
        % create new color, otherwise parcel will not be identifyable
        table_new.table(current_label, 1:3) = table_new.table(current_label, 1:3) + vec; 
        table_new.table(current_label, 5) = table_new.table(current_label, 5) + sum(vec.*[2^0, 2^8, 2^16]);
        label_new(label_alex == current_orig_labels(idx(j))) = table_new.table(current_label, 5);
    end
end
table_new.struct_names = table_new.struct_names';
table_new.numEntries = current_label;
label_new(label_new == 0) = table_new.table(1, 5);

end

function vec = compute_vec_from_label(internal_label)

% This function generates the difference between component color and 
% current parcel's color

% update color. rgb: we have 3 possible binary values, hence this is 8
oct_internal_label = dec2base(abs(internal_label), 8);
vec = 0;
for k = 1:length(oct_internal_label)
    vec(k) = str2num(oct_internal_label(k)); %make a vector form the binary text
end
if(length(vec) == 1) %% make the vector length 3
    vec = [0, 0, vec];
elseif(length(vec) == 2)
    vec = [0, vec];
elseif(length(vec) > 3)
    internal_label
    vec
    error('something went wrong')
end
if(sign(internal_label) == -1)
    vec = -vec;
elseif(sign(internal_label) == 1)
elseif(sign(internal_label) == 0)
else
    internal_label
    error('something went wrong')
end
end

function collision = check_for_collision(vec, new_table, new_label)

% This function checks whether there is a collision between the color of
% current parcel and the colors which are already chosen for previous
% parcels. 

I = ismember(new_table.table(:, 1:3), new_table.table(new_label, 1:3) + vec, 'rows');
collision = sum(I); % collision in all entries
collision = collision + sum(sign(new_table.table(new_label, 1:3) + vec) == -1); % add non zero check
collision = collision + sum((new_table.table(new_label, 1:3) + vec) > 255); % add not above 255 check

end

function off_set_internal_out = check_for_overflow(vec, new_table, new_label, off_set_internal)

% This function generates a new offset when there is a collision

collision = sum((new_table.table(new_label, 1:3) + vec) > 255);% we run out of the color space
if(collision > 0)
    off_set_internal_out = off_set_internal - 256; % go back 10, until we find 'new color space'
else
    off_set_internal_out = off_set_internal;
end

end
