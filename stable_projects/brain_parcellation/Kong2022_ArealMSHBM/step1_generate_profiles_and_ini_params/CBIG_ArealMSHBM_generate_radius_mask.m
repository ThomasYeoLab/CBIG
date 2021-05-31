function CBIG_ArealMSHBM_generate_radius_mask(lh_labels, rh_labels, mesh, radius, out_dir, lh_distance, rh_distance)

% CBIG_ArealMSHBM_generate_radius_mask(lh_labels, rh_labels, mesh, radius, out_dir, lh_distance, rh_distance)
%
% This function will be used to generate radius mask for parcels from a given group-level 
% parcellation. This group-level parcellation (defined by lh_labels and rh_labels) should be 
% the same one used in CBIG_ArealMSHBM_generate_ini_params for initializing the clustering 
% parameters. 
%
% Input:
%
%   - lh_labels rh_labels: (#vertices x 1 vector for both)
%     
%     Labels of a given group-level parcellation. The labels of both lh_labels
%     and rh_labels should starts from 1 to L. The medial wall area should be
%     denoted as 0.
%
%   - mesh: (string)
%     
%     The data surface space. 'fsaverage5/fsaverage6/fsaverage' or 'fs_LR_32k'. 
%
%   - radius: (string)
%
%     The size (in mm) of the radius mask. For example, '30'.
%
%   - out_dir: (string)
%     
%     The spatial radius mask sparse variables <lh_boundary> and <rh_boundary>
%     will be saved in:
%     <out_dir>/spatial_mask/spatial_mask_<mesh>.mat
%
%   - lh_distance, rh_distance (num_verts x num_verts matrices) [optional]
%
%     If user wants to pass in a pre-computed geodesic distance matrices 
%     lh_distance, rh_distance, this script will skip the step of generating
%     distance matrices.  
%
% Examples:
% lh_group_labels = CBIG_read_annotation(['lh.Schaefer2018_400Parcels_17Networks_order.annot']);
% rh_group_labels = CBIG_read_annotation(['rh.Schaefer2018_400Parcels_17Networks_order.annot']);
% lh_labels = lh_group_labels - 1;
% rh_labels = rh_group_labels - 1;
% CBIG_ArealMSHBM_generate_radius_mask(lh_labels, rh_labels, 'fsaverage6', '30', './test_output')
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

radius = str2num(radius);
disp('1. generate geodesic distance matrix for the given mesh ...')
if(nargin == 5)
    disp('==> Left hemi ...')
    lh_distance = generate_geodesic_distance_mat(mesh, 'lh');
    disp('==> Left hemi ... Done')
    disp('==> Right hemi ...')
    rh_distance = generate_geodesic_distance_mat(mesh, 'rh');
    disp('==> Right hemi ... Done')
end
if(nargin == 7)
    disp('Use distance matrices passed in.')
end
disp('2. find the parcels around the central sulcus ...')
[lh_avg_dis, rh_avg_dis] = find_central_sulcus_parcels(lh_distance, rh_distance, mesh, lh_labels, rh_labels);

disp('3. find the radius mask for each parcel ...')
[lh_boundary,rh_boundary] = add_spatial_constraint(lh_distance,rh_distance,lh_labels,rh_labels,mesh,radius);

disp('4. truncate the radius mask so that parcels will not expand across central sulcus ...')
[lh_boundary, rh_boundary] = truncate_spatial_constraint_expand(lh_boundary,...
     rh_boundary, lh_avg_dis, rh_avg_dis, lh_labels, rh_labels, mesh);
lh_boundary = sparse(lh_boundary);
rh_boundary = sparse(rh_boundary);

if(~exist(fullfile(out_dir,'spatial_mask')))
    mkdir(fullfile(out_dir,'spatial_mask'))
end
save(fullfile(out_dir,'spatial_mask',['spatial_mask_' mesh '.mat']),'lh_boundary','rh_boundary');

end

function dist_transform = generate_geodesic_distance_mat(mesh, hemi)

if(~isempty(strfind(mesh, 'fsaverage')))
    MARS_sbjMesh = CBIG_ReadNCAvgMesh(hemi, mesh, 'inflated','cortex');
elseif(~isempty(strfind(mesh, 'fs_LR')))
    MARS_sbjMesh = CBIG_read_fslr_surface(hemi, mesh, 'inflated', 'medialwall.annot');
end
num_vertices = int32(size(MARS_sbjMesh.vertices, 2));
maxNeighbors = int32(size(MARS_sbjMesh.vertexNbors, 1));
dist_transform = zeros(num_vertices, num_vertices);

vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)),...
     int32(size(MARS_sbjMesh.vertices, 2)), int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.vertices));
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);

for i = 1:num_vertices
    
    dist_transform(i, :) = MARS_DT_Boundary(int32([1:1:num_vertices] == i), num_vertices,...
         maxNeighbors, MARS_sbjMesh.vertexNbors, double(vertexDist2Nbors)); %min_heap assumes double
    
end
end

function [lh_avg_dis, rh_avg_dis] = find_central_sulcus_parcels(lh_distance,rh_distance,mesh,lh_labels,rh_labels)
%load mesh
% precentral:25 postcentral:23
if(strcmp(mesh, 'fs_LR_32k'))
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh, 'inflated','medialwall.annot');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh, 'inflated','medialwall.annot');
else
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','aparc.annot');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated','aparc.annot');
end

for l = 1:max(rh_labels)
    lh_avg_dis_precentral(l) = mean(mean(lh_distance(lh_labels == l, lh_avg_mesh.MARS_label == 25)));
    lh_avg_dis_postcentral(l) = mean(mean(lh_distance(lh_labels == l, lh_avg_mesh.MARS_label == 23)));
    
    rh_avg_dis_precentral(l) = mean(mean(rh_distance(rh_labels == l, rh_avg_mesh.MARS_label == 25)));
    rh_avg_dis_postcentral(l) = mean(mean(rh_distance(rh_labels == l, rh_avg_mesh.MARS_label == 23)));
end

% precentral | postcentral
lh_avg_dis = [lh_avg_dis_precentral; lh_avg_dis_postcentral];
rh_avg_dis = [rh_avg_dis_precentral; rh_avg_dis_postcentral];
end

function [lh_boundary,rh_boundary] = add_spatial_constraint(lh_distance,rh_distance,lh_labels,rh_labels,mesh,radius)
if(strcmp(mesh, 'fs_LR_32k'))
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh, 'inflated','medialwall.annot');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh, 'inflated','medialwall.annot');
else
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','aparc.annot');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated','aparc.annot');
end
lh_num_verts = size(lh_avg_mesh.vertices, 2);
rh_num_verts = size(rh_avg_mesh.vertices, 2);

lh_boundary = zeros(lh_num_verts,length(unique(lh_labels))-1);
rh_boundary = zeros(rh_num_verts,length(unique(rh_labels))-1);

%Find boundary vertices for each parcel
for lh_l = unique(lh_labels)'
    if (lh_l ~= 0)
        lh_curr_parcel_verts = lh_labels == lh_l;
        lh_curr_boundary_verts = ConvertObjVec2Boundary(lh_avg_mesh,lh_curr_parcel_verts'); %boundary vertices are 1s
        lh_curr_boundary_dis = lh_distance(lh_curr_boundary_verts,:);
        lh_curr_boundary_dis = lh_curr_boundary_dis <= radius;
        lh_curr_boundary_idx = unique(find(sum(lh_curr_boundary_dis) ~= 0));
        lh_boundary(:,lh_l) = lh_curr_parcel_verts;
        lh_boundary(lh_curr_boundary_idx,lh_l) = 1;   
    end
end
for rh_l = unique(rh_labels)'
    if (rh_l ~= 0)
        rh_curr_parcel_verts = rh_labels == rh_l;
        rh_curr_boundary_verts = ConvertObjVec2Boundary(rh_avg_mesh,rh_curr_parcel_verts'); %boundary vertices are 1s
        rh_curr_boundary_dis = rh_distance(rh_curr_boundary_verts,:);
        rh_curr_boundary_dis = rh_curr_boundary_dis <= radius;
        rh_curr_boundary_idx = unique(find(sum(rh_curr_boundary_dis) ~= 0));
        rh_boundary(:,rh_l) = rh_curr_parcel_verts;
        rh_boundary(rh_curr_boundary_idx,rh_l) = 1;
    end
end
end

function [lh_boundary, rh_boundary] = truncate_spatial_constraint_expand(lh_boundary, rh_boundary,...
         lh_avg_dis, rh_avg_dis, lh_labels, rh_labels, mesh)
%load mesh
% precentral:25 postcentral:23 paracentral:18 
%if we mask out precentral(25), we will further mask superiorfrontal(29),
%caudalmiddlefrontal(4), parsopercularis(19)
%if we mask out postcentral(23), we will further mask superiorparietal(30),
%supramarginal(32)

if(strcmp(mesh, 'fs_LR_32k'))
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh, 'inflated','aparc.annot');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh, 'inflated','aparc.annot');
    lh_aparc_label = lh_avg_mesh.MARS_label;
    rh_aparc_label = rh_avg_mesh.MARS_label;
else
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','aparc.annot');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated','aparc.annot');
    lh_aparc_label = lh_avg_mesh.MARS_label';
    rh_aparc_label = rh_avg_mesh.MARS_label';
end

lh_aparc_label(lh_aparc_label == 4) = 19;
lh_aparc_label(lh_aparc_label == 29) = 19;
lh_aparc_label(lh_aparc_label == 32) = 30;

rh_aparc_label(rh_aparc_label == 4) = 19;
rh_aparc_label(rh_aparc_label == 29) = 19;
rh_aparc_label(rh_aparc_label == 32) = 30;

lh_aparc_verts = zeros(size(lh_aparc_label,1),max(lh_aparc_label));
rh_aparc_verts = zeros(size(rh_aparc_label,1),max(rh_aparc_label));

lh_aparc_verts(sub2ind(size(lh_aparc_verts),[1:1:size(lh_aparc_verts,1)]',lh_aparc_label)) = 1;
rh_aparc_verts(sub2ind(size(rh_aparc_verts),[1:1:size(rh_aparc_verts,1)]',rh_aparc_label)) = 1;


lh_label_verts = zeros(size(lh_labels,1),max(lh_labels));
rh_label_verts = zeros(size(rh_labels,1),max(rh_labels));

lh_verts = 1:1:size(lh_label_verts,1);
rh_verts = 1:1:size(rh_label_verts,1);
lh_label_verts(sub2ind(size(lh_label_verts),lh_verts(lh_labels~=0)',lh_labels(lh_labels~=0))) = 1;
rh_label_verts(sub2ind(size(rh_label_verts),rh_verts(rh_labels~=0)',rh_labels(rh_labels~=0))) = 1;

%compute overlap
lh_overlap = lh_label_verts'*lh_aparc_verts;
rh_overlap = rh_label_verts'*rh_aparc_verts;

lh_overlap = bsxfun(@rdivide,lh_overlap,sum(lh_overlap,2));
rh_overlap = bsxfun(@rdivide,rh_overlap,sum(rh_overlap,2));

%extract parcels lie in para/post/pre-central/insula:
%order as:
%para/post/pre-central/insula/superiorfrontal/caudalmiddlefrontal/parsopercularis/superiofarietal/supramarginal
lh_central = lh_overlap(:,[18,23,25,36,19,30]);
lh_central = bsxfun(@times, lh_central, 1./max(lh_central,[],2));
lh_central(isnan(lh_central)) = 0;
%manually fix
lh_central(lh_central(:,2)==0 & lh_central(:,3)~=0,3) = 1;
lh_central(lh_central(:,2)~=0 & lh_central(:,3)==0,2) = 1;
lh_central_bin = lh_central >= 0.999;

rh_central = rh_overlap(:,[18,23,25,36,19,30]);
rh_central = bsxfun(@times, rh_central, 1./max(rh_central,[],2));
rh_central(isnan(rh_central)) = 0;
rh_central(rh_central(:,2)==0 & rh_central(:,3)~=0,3) = 1;
rh_central(rh_central(:,2)~=0 & rh_central(:,3)==0,2) = 1;
rh_central_bin = rh_central >= 0.999;

lh_central_label = zeros(size(lh_labels,1),1);
rh_central_label = zeros(size(rh_labels,1),1);

for i = 1:6
    lh_parcel_set(i,1) = {find(lh_central_bin(:,i)==1)};
    rh_parcel_set(i,1) = {find(rh_central_bin(:,i)==1)};
    if(i==1 || i==4)
        lh_t = lh_central(lh_central_bin(:,i)==1,:) >0.2;
        rh_t = rh_central(rh_central_bin(:,i)==1,:) >0.2;
        lh_t(:,[5,6]) = 1;
        rh_t(:,[5,6]) = 1;
        % check distance to pre/post-central
        lh_dis = lh_avg_dis(:,lh_central_bin(:,i)==1);
        rh_dis = rh_avg_dis(:,rh_central_bin(:,i)==1);
        lh_pre_or_post = (lh_dis(1,:)-lh_dis(2,:)) < 0; % mask  keep 1
        rh_pre_or_post = (rh_dis(1,:)-rh_dis(2,:)) < 0;
        lh_t(~lh_pre_or_post,2) = 1;
        rh_t(~rh_pre_or_post,2) = 1;
        lh_t(lh_pre_or_post,2) = 0;
        rh_t(rh_pre_or_post,2) = 0;
        lh_t(lh_pre_or_post,3) = 1;
        rh_t(rh_pre_or_post,3) = 1;
        lh_t(~lh_pre_or_post,3) = 0;
        rh_t(~rh_pre_or_post,3) = 0;
        
        lh_parcel_set(i,2) = {lh_t};
        rh_parcel_set(i,2) = {rh_t};

    elseif(i==2)
        lh_t = lh_central(lh_central_bin(:,i)==1,:) >=0.999;
        rh_t = rh_central(rh_central_bin(:,i)==1,:) >=0.999;
        lh_t(:,[1,4,6]) = 1;
        rh_t(:,[1,4,6]) = 1;
        lh_parcel_set(i,2) = {lh_t};
        rh_parcel_set(i,2) = {rh_t};

    elseif(i==3)
        lh_t = lh_central(lh_central_bin(:,i)==1,:) >=0.999;
        rh_t = rh_central(rh_central_bin(:,i)==1,:) >=0.999;
        lh_t(:,[1,4,5]) = 1;
        rh_t(:,[1,4,5]) = 1;
        lh_parcel_set(i,2) = {lh_t};
        rh_parcel_set(i,2) = {rh_t};
        
    elseif(i==5)
        lh_t = ones(size(lh_central_bin));
        rh_t = ones(size(lh_central_bin));
        lh_t(:,2) = 0;
        rh_t(:,2) = 0;
        lh_parcel_set(i,2) = {lh_t};
        rh_parcel_set(i,2) = {rh_t};

    elseif(i==6)
        lh_t = ones(size(lh_central_bin));
        rh_t = ones(size(lh_central_bin));
        lh_t(:,3) = 0;
        rh_t(:,3) = 0;
        lh_parcel_set(i,2) = {lh_t};
        rh_parcel_set(i,2) = {rh_t};  
    end
end

for i = 1:6
    current_mask = lh_parcel_set{i,2};
    current_parcel_set = lh_parcel_set{i,1};
    for j = 1:length(current_parcel_set)
        current_parcel = current_parcel_set(j);
        mask_region = find(current_mask(j,:) == 0);
        for k = mask_region
            if( k==1 ) %mask out paracentral
                lh_boundary(lh_aparc_label == 18, current_parcel) = 0; %paracentral
            elseif( k==2 ) %mask out postcentral
                lh_boundary(lh_aparc_label == 23, current_parcel) = 0; %potcentral
            elseif( k==3 ) %mask out precentral
                lh_boundary(lh_aparc_label == 25, current_parcel) = 0; %precentral
            elseif( k==4 ) %mask out insula
                lh_boundary(lh_aparc_label == 36, current_parcel) = 0; %insula
            elseif( k==5 ) %mask out superiorfrontal/caudalmiddlefrontal/parsopercularis
                lh_boundary(lh_aparc_label == 19, current_parcel) = 0; 
            elseif( k==6 ) %mask out superiofarietal/supramarginal
                lh_boundary(lh_aparc_label == 30, current_parcel) = 0;
            end
        end
    end
end
            
for i = 1:6
    current_mask = rh_parcel_set{i,2};
    current_parcel_set = rh_parcel_set{i,1};
    for j = 1:length(current_parcel_set)
        current_parcel = current_parcel_set(j);
        mask_region = find(current_mask(j,:) == 0);
        for k = mask_region
            if( k==1 ) %mask out paracentral
                rh_boundary(rh_aparc_label == 18, current_parcel) = 0; %paracentral
            elseif( k==2 ) %mask out postcentral
                rh_boundary(rh_aparc_label == 23, current_parcel) = 0; %potcentral
            elseif( k==3 ) %mask out precentral
                rh_boundary(rh_aparc_label == 25, current_parcel) = 0; %precentral
            elseif( k==4 ) %mask out insula
                rh_boundary(rh_aparc_label == 36, current_parcel) = 0; %insula
            elseif( k==5 ) %mask out superiorfrontal/caudalmiddlefrontal/parsopercularis
                rh_boundary(rh_aparc_label == 19, current_parcel) = 0; 
            elseif( k==6 ) %mask out superiofarietal/supramarginal
                rh_boundary(rh_aparc_label == 30, current_parcel) = 0;
            end
        end
    end
end

end                