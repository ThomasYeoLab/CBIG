function CBIG_SPGrad_generate_gradient_matrix(mesh, medial_mask, downsample, output_dir)

% CBIG_SPGrad_generate_gradient_matrix(mesh, medial_mask, downsample, output_dir)
%
% This script reads in the resting-state function connectivity gradient map stored in 
% <output_dir>/gradients_edge_density.dtseries.nii, downsamples the gradient map, and
% computes the pairwise distance matrix bwteen vertices. The distance between vertcies
% is defined as the shortest path between them, where the edge weight is defined by
% the gradient value.
%
% Input:
%     - mesh:
%       resolution of surface mesh, e.g. 'fsaverage6', 'fs_LR_32k'
%
%     - medial_mask: (#num_vertices x 1 binary vector or 'NONE')
%       the medial area mask. Set it as 'NONE' unless the medial wall area is defined differently
%       as fsaverage medial wall area from FreeSurfer or fs_LR_32k medial wall area from HCP. 
%       If the data has a different medial area mask, please pass in a #num_vertices x 1 binary
%       vector. The medial area should be denoted as 1, others should be denoted as 0.
%
%     - downsample: (string)
%       set downsample to be '0' to not perform downsampling on the gradient map. Set downsample
%       to be a positive number, the gradient map will be downsampled to #vertices/(2*<downsample>) 
%       For fsaverage6, we suggest user set: downsample = '3.2'.
%       For fs_LR_32k, we suggest user set: downsample = '3.2'.
%
%     - output_dir:
%       output directory to save the results. The estimated gradient distance matrix will be 
%       saved as:
%       <output_dir>/lh_gradient_distance_matrix.npy
%       <output_dir>/rh_gradient_distance_matrix.npy
%
% Example:
%   CBIG_SPGrad_generate_gradient_matrix('fsaverage6', 'NONE', '3.2', out_dir);
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'),...
 '/external_packages/matlab/non_default_packages/cifti-matlab-WashU-gradient')));

downsample = str2num(downsample);

if(strcmp(mesh,'fs_LR_32k'))
    lh_surf = CBIG_read_fslr_surface('lh','fs_LR_32k','sphere','medialwall.annot');
    rh_surf = CBIG_read_fslr_surface('rh','fs_LR_32k','sphere','medialwall.annot');
    if(strcmp(medial_mask,'NONE'))
        medial_mask = [lh_surf.MARS_label == 1;rh_surf.MARS_label == 1];
    end
elseif(~isempty(strfind(mesh,'fsaverage')))
    lh_surf = CBIG_ReadNCAvgMesh('lh', mesh, 'sphere', 'cortex');
    rh_surf = CBIG_ReadNCAvgMesh('rh', mesh, 'sphere', 'cortex');
    if(strcmp(medial_mask,'NONE'))
        medial_mask = [lh_surf.MARS_label==1 rh_surf.MARS_label==1]';
    end
end
input_edge = ft_read_cifti_mod([output_dir '/gradients_edge_density.dtseries.nii']);
input_grad_edge(~medial_mask,1) = input_edge.data;

disp('-----------------------------------------------------------------')
disp('###########################')
disp('#Downsample gradient edges#')
disp('###########################')

if(downsample == 0)
    disp('Do not perform downsampling ...')
    lh_grad_edge = input_grad_edge(1:size(lh_surf.vertices,2));
    rh_grad_edge = input_grad_edge(size(lh_surf.vertices,2)+1:end);
    clear input_grad_edge
else
    disp('Perform downsampling ...')

    downsample = round(size(input_grad_edge,1)/(2*downsample));

    outputtmpdir = [output_dir '/tmp'];

    if (~exist(outputtmpdir))
        mkdir(outputtmpdir);
    end

    if(exist([outputtmpdir '/downsampled.lh.sphere']))
        disp('---------------------> Downsample sphere already exist. Skip ...')
    else
        system(['wb_command -surface-create-sphere ', num2str(downsample), ' ',...
         outputtmpdir, '/sphere.downsampled.R.surf.gii']);
        system(['wb_command -surface-create-sphere ', num2str(downsample), ' ',...
         outputtmpdir, '/sphere.downsampled.L.surf.gii']);
        system(['wb_command -set-structure', ' ', outputtmpdir, '/sphere.downsampled.R.surf.gii CORTEX_RIGHT']);
        system(['wb_command -set-structure', ' ', outputtmpdir, '/sphere.downsampled.L.surf.gii CORTEX_LEFT']);
        surf_gii = gifti([outputtmpdir '/sphere.downsampled.R.surf.gii']);
        write_surf(surf_gii.vertices, surf_gii.faces-1, [outputtmpdir '/downsampled.rh.sphere']);
        surf_gii = gifti([outputtmpdir '/sphere.downsampled.L.surf.gii']);
        write_surf(surf_gii.vertices, surf_gii.faces-1, [outputtmpdir '/downsampled.lh.sphere']);
    end

    lh_down_surf = CBIG_SPGrad_read_surf_mesh([outputtmpdir, '/downsampled.lh.sphere']);
    rh_down_surf = CBIG_SPGrad_read_surf_mesh([outputtmpdir, '/downsampled.rh.sphere']);

    
    [lh_grad_edge, ~, ~] = MARS_linearInterpolate(lh_down_surf.vertices, lh_surf, ...
    input_grad_edge(1:size(lh_surf.vertices,2))');
    [rh_grad_edge, ~, ~] = MARS_linearInterpolate(rh_down_surf.vertices, rh_surf, ...
    input_grad_edge(size(lh_surf.vertices,2)+1:end)');
    clear input_grad_edge
    
    lh_surf = lh_down_surf;
    rh_surf = rh_down_surf;
        
end

disp('Next step ...')
disp('#########################')
disp('#Compute distance matrix#')
disp('#########################')

disp('Construct graph ...')
tic;
lh_G = CBIG_SPGrad_create_graph(lh_surf, lh_grad_edge);
rh_G = CBIG_SPGrad_create_graph(rh_surf, rh_grad_edge);
toc;

disp('Compute pairwise shortest distance ...')
tic;
lh_dist = distances(lh_G);
rh_dist = distances(rh_G);
toc;
lh_dist = single(lh_dist/max(max(abs(lh_dist))));
rh_dist = single(rh_dist/max(max(abs(rh_dist))));
    
% save as NPY format
disp('Save as NPY format ...')            
shape = size(lh_dist);
dataType = class(lh_dist);
header = constructNPYheader(dataType, shape);
fid = fopen([output_dir, '/lh_gradient_distance_matrix.npy'], 'w');
fwrite(fid, header, 'uint8');
fwrite(fid, lh_dist, dataType);
fclose(fid);
clear lh_dist

fid = fopen([output_dir, '/rh_gradient_distance_matrix.npy'], 'w');
fwrite(fid, header, 'uint8');
fwrite(fid, rh_dist, dataType);
fclose(fid);
clear rh_dist

disp('########## Done!') 

rmpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), ...
'/external_packages/matlab/non_default_packages/cifti-matlab-WashU-gradient')));



    
    
