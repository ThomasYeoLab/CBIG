function [lh_rot, rh_rot] = CBIG_spin_test_rotation_surf(mesh, num_rotations, seed)
% [lh_rot, rh_rot] = CBIG_spin_test_rotation_surf(mesh, num_rotations, seed)
%
% This function is used to compute the new vertex points of the left and right
% hemisphere after rotation, relative to the original vertex points in the mesh.
% To compute this, the function will obtain the vertices of the surface mesh,
% and then obtain the rotation matrix. Finally, linear interpolation will be
% performed using the vertices of the surface mesh and rotation matrix as inputs.
%
% INPUT:
%       - mesh:
%         'fs_LR_32k','fsaverage5','fsaverage6', or 'fsaverage'.
%
%       - num_rotations
%         A scalar value used to indicate the desired number of rotations.
%         The default value is 1000.
%
%       - seed
%         A scalar value used to define the random seed for the purpose of
%         producing random rotation.
%         The default value is 1.
%
% OUTPUT:
%       - lh_rot, rh_rot:
%         A num_rotation x #vertices matrix, the vertex points after
%         rotation, relative to the original vertex points in the mesh
%
% Example:
% mesh = 'fsaverage6';
% num_rotations = 5000;
% seed = 2;
% [lh_rot, rh_rot] = CBIG_spin_test_rotation_surf(mesh, num_rotations, seed);
% lh_labels_r = lh_labels(lh_rot(n,:));
% rh_labels_r = rh_labels(rh_rot(n,:));
%Written by Trevor Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Defining the Default Values

    if nargin == 1 % If no input for num_rotations, set default num_rotations as 1000
        num_rotations = 1000;
    end

    if nargin <= 2 % If no input for seed, set default seed as 1
        seed = 1;
    end

%% if fsLR space

    if strcmp(mesh,'fs_LR_32k') == 1
        % Reading in surface mesh for left hemisphere
        lh_avg_mesh = CBIG_read_fslr_surface('lh',mesh,'sphere','medialwall.annot');
        % Reading in surface mesh for right hemisphere
        rh_avg_mesh = CBIG_read_fslr_surface('rh',mesh,'sphere','medialwall.annot');
        verticesl = lh_avg_mesh.vertices; % Reading in vertices of mesh for left hemisphere
        verticesr = rh_avg_mesh.vertices; % Reading in vertices of mesh for right hemisphere
                        
        % Executing CBIG_uniform_rand_rotation function to obtain
        % 3x3xnum_rotation matrix          
        rotation_mat = CBIG_uniform_rand_rotation(num_rotations,seed);

        for i = 1:size(rotation_mat,3)
            verticesl_rot = rotation_mat(:,:,i)*verticesl; % Obtain vertex points on left hemisphere mesh post-rotation
            verticesr_rot = rotation_mat(:,:,i)*verticesr; % Obtain vertex points on right hemisphere mesh post-rotation
            
            % Execution of linear interpolation function to obtain interpolated
            % value(s) of data at vertex points on left hemisphere mesh
            lh_rot(i,:) = MARS_NNInterpolate_kdTree(verticesl_rot, lh_avg_mesh, 1:1:length(verticesl));
            
            % Execution of linear interpolation function to obtain interpolated
            % value(s) of data at vertex points on right hemisphere mesh
            rh_rot(i,:) = MARS_NNInterpolate_kdTree(verticesr_rot, rh_avg_mesh, 1:1:length(verticesr));    
        end
    end

%% if fsaverage space

    if contains(mesh,'fsaverage') == 1
        lh_avg_mesh = CBIG_ReadNCAvgMesh('lh',mesh,'sphere','cortex'); % Reading in surface mesh for left hemisphere
        rh_avg_mesh = CBIG_ReadNCAvgMesh('rh',mesh,'sphere','cortex'); % Reading in surface mesh for left hemisphere
        verticesl = lh_avg_mesh.vertices; % Reading in vertices of mesh for left hemisphere
        verticesr = rh_avg_mesh.vertices; % Reading in vertices of mesh for right hemisphere
             
        % Executing CBIG_uniform_rand_rotation function to obtain
        % 3x3xnum_rotation matrix  
        rotation_mat = CBIG_uniform_rand_rotation(num_rotations,seed);

        for i = 1:size(rotation_mat,3)
            verticesl_rot = rotation_mat(:,:,i)*verticesl; % Reading in vertices of mesh for left hemisphere
            verticesr_rot = rotation_mat(:,:,i)*verticesr; % Reading in vertices of mesh for right hemisphere
            
            % Execution of linear interpolation function to obtain interpolated
            % value(s) of data at vertex points on left hemisphere mesh
            lh_rot(i,:) = MARS_NNInterpolate_kdTree(verticesl_rot, lh_avg_mesh, 1:1:length(verticesl));
            
            % Execution of linear interpolation function to obtain interpolated
            % value(s) of data at vertex points on right hemisphere mesh
            rh_rot(i,:) = MARS_NNInterpolate_kdTree(verticesr_rot, rh_avg_mesh, 1:1:length(verticesr));      
        end        
    end
end