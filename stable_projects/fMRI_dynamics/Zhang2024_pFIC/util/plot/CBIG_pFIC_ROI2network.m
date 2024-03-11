function network_assignment = CBIG_pFIC_ROI2network(input_roi)

% network_assignment = CBIG_pFIC_ROI2network(input_roi)
% This function first upsamples ROI-level values to vertex, then assign
% vertices to yeo 7 networks with the maximal overlapping. This function is mainly used to generate the
% boxplot depicting the network pattern of the input_roi.
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Note that this function only converts from Deiskan parcellation to yeo 7 network
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% Input:
%   - input_roi: a 72-by-1 vector arranged according to Desikan
%   parcellation. (See ../general/Deiskan_72.txt for more details)
% Output:
%   - network_assignment: a 14114-by-2 matrix. The first column is the network index 
% from 1 to 7. The second column is the value of vertices corresponding to
% the network.
%
% Example:
% network_assignment = CBIG_pFIC_ROI2network(EI_contrast);
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%%  Desikan parcellation with 72 cortical ROI, See 'desikan_72.txt'
% settiing 4 ROIs corresponding to medial wall to 0 (in case they are not already 0)
input_roi([1 5 37 41]) = 0;

%% load network definition
network_def = load('yeo7network_51_roi_label.mat');
lh_label_network = network_def.lh_labels;
rh_label_network = network_def.rh_labels;
lh_label_network_roi_name = network_def.lh_labels_name;
rh_label_network_roi_name = network_def.rh_labels_name;

% initialize vertex in fsaverage5 space
lh_label_network_fs5 = zeros(10242,1);
rh_label_network_fs5 = zeros(10242,1);

lh_label_network_fs5(lh_label_network(:,1)) = lh_label_network(:,2);
rh_label_network_fs5(rh_label_network(:,1)) = rh_label_network(:,2);
% add an offset to differentiate rh from lh labels
rh_label_network_fs5(rh_label_network_fs5>0) = rh_label_network_fs5(rh_label_network_fs5>0) + 26; 

%% read Desikan labels
lh_mesh_fs5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'aparc.annot');
rh_mesh_fs5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'aparc.annot');

lh_label_desikan = lh_mesh_fs5.MARS_label';
rh_label_desikan = rh_mesh_fs5.MARS_label';

for i = 1:length(input_roi)/2
       if input_roi(i) == 0
           lh_label_desikan(lh_label_desikan == i) = 0;
       end
end
for i = length(input_roi)/2:length(input_roi)
       if input_roi(i) == 0
           rh_label_desikan((rh_label_desikan + 36) == i) = 0;
       end
end

%% upsample the ROI-level parameter values to vertex-level
lh_param = input_roi(1:length(input_roi)/2);
rh_param = input_roi(length(input_roi)/2+1:length(input_roi));
lh_param_vertex = zeros(10242,1);
rh_param_vertex = zeros(10242,1);

for i = 1:10242
    if lh_label_desikan(i) == 0
        lh_param_vertex(i) = NaN;
    else
        lh_param_vertex(i) = lh_param(lh_label_desikan(i));
    end
end
for i = 1:10242
    if rh_label_desikan(i) == 0
        rh_param_vertex(i) = NaN;
    else
        rh_param_vertex(i) = rh_param(rh_label_desikan(i));
    end
end


%% group vertices according to Yeo 51 ROIs (components)
for i = 1:26 
    lh_param_vertex_component = lh_param_vertex((lh_label_network_fs5 == i));
    lh_param_vertex_component(isnan(lh_param_vertex_component)) = []; 
    lh_rh_param_vertex_component{i} = lh_param_vertex_component;
end
for i = 27:51 
    rh_param_vertex_component = rh_param_vertex((rh_label_network_fs5 == i));
    rh_param_vertex_component(isnan(rh_param_vertex_component)) = [];  
    lh_rh_param_vertex_component{i} = rh_param_vertex_component;
end

%% assign 51 components to 7 networks
% 7 network names
network_name = {'SomMot'; 'Vis'; 'DorsAttn'; 'Sal/VentAttn'; 'Limbic'; 'Control'; 'Default'};
label_51_to_7 = zeros(51,1);
label_51_to_7(1) = 2;
label_51_to_7(2) = 1;
label_51_to_7(3:5) = 3;
label_51_to_7(6:10) = 4;
label_51_to_7(11:12) = 5;
label_51_to_7(13:21) = 6;
label_51_to_7(22:26) = 7;
label_51_to_7(27) = 2;
label_51_to_7(28) = 1;
label_51_to_7(29:31) = 3;
label_51_to_7(32:37) = 4;
label_51_to_7(38:39) = 5;
label_51_to_7(40:46) = 6;
label_51_to_7(47:51) = 7;

network_param = [];
network_index = [];

for i = 1:7
    component_index = find(label_51_to_7 == i);
    for j = 1:length(component_index)
        network_param = [network_param; lh_rh_param_vertex_component{component_index(j)}];
        % from 1 to 7
        network_index = [network_index; i*ones(length(lh_rh_param_vertex_component{component_index(j)}), 1)]; 
    end
end

network_assignment(:, 1) = network_index;
network_assignment(:, 2) = network_param;

end
