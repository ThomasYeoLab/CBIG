function CBIG_Yeo2011_GrowBoundaries(num_network, outdir, save_int_file) 

% CBIG_Yeo2011_GrowBoundaries(num_network, outdir, save_int_file)
%
% This function takes the Yeo2011 Split Components and grow the boundaries.
% The Yeo2011 Split Components Label has its boundaries being eroded. This
% code aims to reverse the erosion by growing back the boundaries and tries
% to match the Yeo2011 7 or 17 Networks labels as much as possible.This
% single function will grow back both the rh and lh hemisphere boundaries.
% Note that this function assumes that the Yeo2011 files are stored in 
% [CBIG_CODE_DIR '/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/
% 1000subjects_reference/Yeo_JNeurophysiol11_SplitLabels/fsaverage5/label/'];
%
% Inputs:
%   - num_network
%     A binary scalar or string of 7 or 17.The user will use this input
%     variable to define if they want to grow back the 7 or 17 Networks
%     Split Components label.
%
%   - outdir
%     A string. This string defines the output directory where the user
%     wants to save the growed boundaries in. The final output files will
%     be saved in [outdir '/rh.growed_final.annot'] for rh hemisphere and
%     [outdir '/lh.growed_final.annot'] for lh hemisphere.
%
%   - save_int_file (optional)
%     A binary scalar of 1 or 0. If the user defined this input variable as
%     1, then intermediate files will be saved. If 0 is defined, then no
%     intermediate files will be saved. Saving intermediate files may help
%     in troubleshooting if necessary.If input variable is empty, then it is
%     assumed that save_int_file = 0. As there are 3 steps to this growing
%     process, for step 1, '/*.growed_step1.annot' will be saved, where *
%     refers to 'lh' or 'rh'. Step 2 involves interations and hence the
%     following files will be saved '/*h.growed_' iter_num '.annot'. Step 3
%     is the final step and hence the final growed file will be saved (as
%     defined in the 'outdir' description.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Yan Rui Tan

%% Setup
if ~ischar(num_network)
    num_network = num2str(num_network);
end

if ~exist('save_int_file')
    save_int_file = 0;
else
    if ischar(save_int_file) == 1
        save_int_file = str2num(save_int_file);
    end
end

mkdir(outdir)

%% Load Split Component Labels
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

annot_path = fullfile(CBIG_CODE_DIR,'stable_projects','brain_parcellation', ...
    'Yeo2011_fcMRI_clustering','1000subjects_reference', ...
    'Yeo_JNeurophysiol11_SplitLabels', 'fsaverage5', 'label');

lh_annot = fullfile(annot_path, ...
    ['lh.Yeo2011_' num_network 'Networks_N1000.split_components.annot']);
rh_annot = fullfile(annot_path, ...
    ['rh.Yeo2011_' num_network 'Networks_N1000.split_components.annot']);


[lh_vertices, lh_label, lh_colortable] = read_annotation(lh_annot);
[rh_vertices, rh_label, rh_colortable] = read_annotation(rh_annot);
 
%% Remove medial wall from neighbors matrix
lh_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
rh_mesh = CBIG_ReadNCAvgMesh('rh','fsaverage5','inflated','cortex');

lh_medial_label = find(lh_mesh.MARS_label==1); 
rh_medial_label = find(rh_mesh.MARS_label==1);

lh_label(lh_medial_label) = -1;
rh_label(rh_medial_label) = -1;

lh_neighbor = lh_mesh.vertexNbors;

lh_neighbor_reshape = reshape(lh_neighbor, size(lh_neighbor,2) * 6, 1);
for i = 1:length(lh_medial_label)
    curr_med_label = lh_medial_label(i);
    a = find(lh_neighbor_reshape == curr_med_label);
    lh_neighbor_reshape(a) = 0;
end

lh_neighbor = reshape(lh_neighbor_reshape, 6, size(lh_neighbor,2));

rh_neighbor = rh_mesh.vertexNbors;

rh_neighbor_reshape = reshape(rh_neighbor, size(rh_neighbor,2) * 6, 1);
for i = 1:length(rh_medial_label)
    curr_med_label = rh_medial_label(i);
    a = find(rh_neighbor_reshape == curr_med_label);
    rh_neighbor_reshape(a) = 0;
end

rh_neighbor = reshape(rh_neighbor_reshape, 6, size(rh_neighbor,2));

%% Change label in 17 Networks Split Components
% There are 2 vertices that were wrongly labelled in 17 Networks Split
% Components. This can be seen by comparing the Split Component Version and
% the original 17 Networks Version. By right, the 2 vertices should be labelled
% as the same network. However, they are labelled differently in
% the split components and the original 17Networks version. So to correct for it,
% the split components label for that 2 vertices is manually changed.

if str2num(num_network) == 17
    lh_label(9247) = 16399300;
    lh_label (9246) = 16399300;
end

%% Perform Step 1
% In between 2 parcels are separated by black boundaries that are at least
% 2 vertices thick. So the first step involves growing each parcel by 1
% parcel from its boundary.

% Starting with Lh
tmp_unique = unique(lh_label);
tmp_unique(tmp_unique == -1) = [];
tmp_unique(tmp_unique == 65793) = [];

for i = 1:length(tmp_unique)
    curr_unique = tmp_unique(i);
    curr_find_all = find(lh_label == curr_unique);
    curr_all_neigh = lh_neighbor(:,curr_find_all);
    curr_all_neigh = reshape(curr_all_neigh, numel(curr_all_neigh),1);
    curr_all_neigh = curr_all_neigh(curr_all_neigh ~= 0);
    
    neigh_labels = lh_label(curr_all_neigh);
    curr_all_neigh = curr_all_neigh(neigh_labels ~= curr_unique);
    lh_label(curr_all_neigh) = curr_unique;
end

% Then with RH

tmp_unique = unique(rh_label);
tmp_unique(tmp_unique == -1) = [];
tmp_unique(tmp_unique == 65793) = [];

for i = 1:length(tmp_unique)
    curr_unique = tmp_unique(i);
    curr_find_all = find(rh_label == curr_unique);
    curr_all_neigh = rh_neighbor(:,curr_find_all);
    curr_all_neigh = reshape(curr_all_neigh, numel(curr_all_neigh),1);
    curr_all_neigh = curr_all_neigh(curr_all_neigh ~= 0);
    
    neigh_labels = rh_label(curr_all_neigh);
    curr_all_neigh = curr_all_neigh(neigh_labels ~= curr_unique);
    check = length(unique(neigh_labels));
    rh_label(curr_all_neigh) = curr_unique;
end


%% Write annotation files for step 1
if save_int_file == 1
    write_annotation(fullfile(outdir, 'lh.growed_step1.annot'), lh_vertices, ...
        lh_label, lh_colortable);
    write_annotation(fullfile(outdir, 'rh.growed_step1.annot'), rh_vertices, ...
        rh_label, rh_colortable);
end

%% Preparing for Step 2
% As we want to grow the boundaries such that the final growed version is
% very similar to the Yeo2011 7 or 17 Networks, we need to use the Yeo2011
% 7 or 17 Networks file to grow the boundaries.

lh_true_annot_path = fullfile(annot_path,...
    ['lh.Yeo2011_' num_network 'Networks_N1000.annot']);
rh_true_annot_path = fullfile(annot_path, ...
    ['rh.Yeo2011_' num_network 'Networks_N1000.annot']);

[lh_true_vertices, lh_true_label, lh_true_colortable] = read_annotation(lh_true_annot_path);
[rh_true_vertices, rh_true_label, rh_true_colortable] = read_annotation(rh_true_annot_path);

% Because the Split Component labels and the 7/17 Networks labels are
% different. In order to use the 7/17 Networks labels to help us grow
% boundaries, we need to match the labels between the 2 types of labels.
lh_matched_label = match_network_labels (lh_true_annot_path, lh_label, lh_colortable);
rh_matched_label = match_network_labels (rh_true_annot_path, rh_label, rh_colortable);


%% Step 2
% Because the Yeo2011 7 or 17 Networks have some small isolated network
% parcels (ie group of very few vertices assigned to some networks that is
% intersperse among other big networks). So when we grow, we do not want to
% grow these small network parcels back. So how this works is we first find
% the black vertices in the split component version.
% We go in order of the vertices that have the most number of neighbors
% that are labelled. For each vertex, we check its network label 
% in the 7/17 Network version. For example, we find that it is Network A.
% Then we check in the split component parcellation the component labels of
% its neighbors. From the component labels, if we find that we have at
% least one that come from Network A, then we will assign the vertex to be
% the label from the components from Network A. If the vertex neighbors
% have more than 1 components from Network A, then we will see the number
% of neighbors from each component and assign the vertex to be the
% component label with the more neighbours.

lh_bord_vert = find(lh_label == 65793);
rh_bord_vert = find(rh_label == 65793);
check_lh = 1;
check_rh = 1;

count_lh = 0;
count_rh = 0;

% Start lh iteration
while check_lh ~= 0
    
    [lh_unique_color, lh_all_colors, ~, lh_nb_unique] = ...
        Summarise(lh_bord_vert, lh_neighbor, lh_label);    
    max_lh_nb_unique = max(lh_nb_unique);   
    num_lh_bord_init = length(lh_bord_vert);    
    [lh_label] = Grow_Step2(lh_unique_color, lh_bord_vert, lh_nb_unique,...
        max_lh_nb_unique, lh_label, lh_all_colors, lh_neighbor,...
        lh_true_label, lh_matched_label, lh_colortable.table(:,5));    
    lh_bord_vert = find(lh_label == 65793);        
    check_lh = num_lh_bord_init - length(lh_bord_vert);        
    count_lh = count_lh + 1;   
    if save_int_file == 1
        write_annotation(fullfile(outdir, ['lh.growed_' num2str(count_lh) '.annot']),...
        lh_vertices, lh_label, lh_colortable);
    end
    
end

% Start rh iterations
while check_rh ~= 0
    
    num_rh_bord_init = length(rh_bord_vert);
    [rh_unique_color, rh_all_colors, ~, rh_nb_unique] = ...
        Summarise(rh_bord_vert, rh_neighbor, rh_label);    
    max_rh_nb_unique = max(rh_nb_unique);    
    [rh_label] = Grow_Step2(rh_unique_color, rh_bord_vert, rh_nb_unique, ...
        max_rh_nb_unique, rh_label, rh_all_colors, rh_neighbor, rh_true_label,...
        rh_matched_label, rh_colortable.table(:,5));    
    rh_bord_vert = find(rh_label == 65793);    
    check_rh = num_rh_bord_init - length(rh_bord_vert);       
    count_rh = count_rh + 1;
    if save_int_file == 1
        write_annotation(fullfile(outdir, ['rh.growed_' num2str(count_rh) '.annot']),...
        rh_vertices, rh_label, rh_colortable);
    end
     
end


%% Step 3
% After step 2, what should remain are the black vertices that originally
% belonged to small Network Parcels. Now, we grow these black vertices,
% also starting from those vertices that have more neighbours that are
% labelled.To assign the component label of that vertex we now check the
% component labels of the neghbors. Whichever component that is assigned to
% the most numbr of neihbors, the vertex will similarly be assigned to that
% component label.

lh_bord_vert = find(lh_label == 65793);
rh_bord_vert = find(rh_label == 65793);
check_lh = 1;
check_rh = 1;

% Start lh iteration
while check_lh ~=0
    
    [~, lh_all_colors, ~, lh_nb_unique] = ...
        Summarise(lh_bord_vert, lh_neighbor, lh_label);   
    max_lh_nb_unique = max(lh_nb_unique);            
    [lh_label] = Grow_Step3(lh_bord_vert, lh_nb_unique, max_lh_nb_unique,...
        lh_label, lh_all_colors);   
    lh_bord_vert = find(lh_label == 65793);    
    if isempty(lh_bord_vert)
        check_lh = 0;
    end
    
end

% Start rh iterations
while check_rh ~=0
    
    [~, rh_all_colors, ~, rh_nb_unique] = ...
        Summarise(rh_bord_vert, rh_neighbor, rh_label);   
    max_rh_nb_unique = max(rh_nb_unique);             
    [rh_label] = Grow_Step3(rh_bord_vert, rh_nb_unique, max_rh_nb_unique,...
        rh_label, rh_all_colors);   
    rh_bord_vert = find(rh_label == 65793);   
    if isempty(rh_bord_vert)
        check_rh = 0;
    end
    
end

% Convert medial wall labels back to black
lh_label(lh_medial_label) = 65793;
rh_label(rh_medial_label) = 65793;

% Save final output file
write_annotation(fullfile(outdir,...
    ['lh.Yeo2011_' num_network 'Networks_growed.annot']), lh_vertices,...
    lh_label, lh_colortable);
write_annotation(fullfile(outdir ,...
    ['rh.Yeo2011_' num_network 'Networks_growed.annot']), rh_vertices,...
    rh_label, rh_colortable);


function assigned_label = match_network_labels (true_annot_path, label, colortable)

[~, true_label, true_colortable] = read_annotation(true_annot_path);

true_col = true_colortable.table(:,5);
col = colortable.table(:,5);

true_col(true_col == 65793) = [];
col (col == 65793) = [];

for w = 1:length(true_col)
    curr_col = true_col(w);
    true_vertex{w} = find(true_label == curr_col);
end

assigned_label(1,:) = 65793;
for w = 1:length(col)
    curr_col = col(w);
    vertex = find(label == curr_col);
    count = 0;
    for j = 1:length(true_col)
        curr_true_vertex = true_vertex{j};
        match = intersect(vertex, curr_true_vertex);
        num_match = length(match);
        if num_match > count
            count = num_match;
            assigned_label(w+1,:) = true_col(j);
        end
    end
end

function [label] = Grow_Step2(unique_color, bord_vert, nb_unique, max_nb_unique,...
        label, all_colors, neighbor, true_label, assigned_label, test_label)

check1 = 1;
while check1 ~= 0
    tmp_find1 = find(nb_unique == max_nb_unique);

    for p = 1:length(tmp_find1)

        curr_value1 = tmp_find1(p);
        curr_color1 = all_colors(:,curr_value1);
        curr_color1(curr_color1 == 65793 ) = [];
        curr_idx1 = bord_vert(curr_value1);

        curr_neighbor = neighbor(:,curr_idx1);
        curr_neighbor = curr_neighbor(curr_neighbor ~= 0);
        curr_true_neighbor_label = true_label(curr_neighbor);
        count_tmp = 0;
        
        % true label is the label in the original 17 Networks of all
        % vertex.So check_label is the corresponding component labels in
        % that network.
        check_label = test_label(assigned_label == true_label(curr_idx1));
        % curr_colo1 stores all the component labels of the neighbors of curr_idx1
        % We check if there are any overlapping component labels between
        % check_label and curr_color1.
        tmp_true_label = intersect(check_label, curr_color1);
        if isempty(tmp_true_label)
            count_tmp = 1;
        end
        
        % If there is overlap, we assign the component labels
        if count_tmp ~= 1
            if length(tmp_true_label) > 1
                count1 = 0;
                for k = 1: length(tmp_true_label)
                    curr_tmp_label = tmp_true_label(k);
                    tmp2 = length(find(curr_color1 == curr_tmp_label));
                    if tmp2 > count1
                        most_color = curr_tmp_label;
                        count1 = tmp2;
                    end
                end
            else
                most_color = tmp_true_label;
            end

            label(curr_idx1) = most_color;
            check1 = 0;
            
        end    

    end

    max_nb_unique = max_nb_unique - 1;
    if max_nb_unique == 0
        check1 = 0;
    end
end


function [label] = Grow_Step3( bord_vert, nb_unique, max_nb_unique, label,...
    all_colors)

tmp_find = find(nb_unique == max_nb_unique);

for m = 1:length(tmp_find)

    curr_value = tmp_find(m);
    curr_color = all_colors(:,curr_value);
    curr_color(curr_color == 65793 ) = [];
    curr_idx = bord_vert(curr_value);
    curr_vertex = mode(curr_color);
    label(curr_idx) = curr_vertex;
end

function [unique_colour,all_colors, num_per_color, non_black_unique] = Summarise(bord_vert, neighbor, label)

unique_colour = zeros(length(bord_vert),6);
num_per_color = zeros(length(bord_vert),6);

for j = 1:length(bord_vert)
    
    tmp_bord_vert = bord_vert(j);
    tmp_neigh = neighbor(:,tmp_bord_vert);
    for i = 1:length(tmp_neigh)
        curr_tmp_neigh = tmp_neigh(i);
        if curr_tmp_neigh > 0
            tmp_label(i,1) = label(curr_tmp_neigh);
        else
            tmp_label(i,1) = nan;
        end
    end
    all_colors(:,j) = tmp_label;
    tmp = unique(tmp_label(~isnan(tmp_label)));

    num_unique = length(tmp);
    unique_colour(j,1:num_unique) = tmp';
    
    for k = 1:num_unique
        total_unique = sum(tmp_label == tmp(k));
        num_per_color(j,k) = total_unique;
    end
    non_black_unique(j,1) = sum(tmp_label(~isnan(tmp_label)) ~= 65793);
end

