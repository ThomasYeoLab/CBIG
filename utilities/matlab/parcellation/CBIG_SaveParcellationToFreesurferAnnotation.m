function CBIG_SaveParcellationToFreesurferAnnotation(matlab_parcellation_file, lh_output_file, rh_output_file)

% CBIG_SaveParcellationToFreesurferAnnotation(matlab_parcellation_file, lh_output_file, rh_output_file)
%
% Wrapper function to write colors of parcellation into FreeSurfer
% annotation file.
% 
% "matlab_parcellation_file" should have "lh_labels" and "rh_labels"
% variables inside.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md




lh_avg_mesh5 = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'sphere', 'cortex');
rh_avg_mesh5 = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'sphere', 'cortex');
lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'sphere', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'sphere', 'cortex');
load freesurfer_color.mat

load(matlab_parcellation_file);
lh_labels7 = MARS_NNInterpolate(lh_avg_mesh.vertices, lh_avg_mesh5, lh_labels');
rh_labels7 = MARS_NNInterpolate(rh_avg_mesh.vertices, rh_avg_mesh5, rh_labels');

num_clusters = length(unique([lh_labels; rh_labels]));

color = freesurfer_color(1:num_clusters, :)*255;
color(1, :) = 1; % freesurfer does not like 0 0 0

CBIG_WriteParcellationToAnnotation(lh_labels7, lh_output_file, color);
CBIG_WriteParcellationToAnnotation(rh_labels7, rh_output_file, color);

exit
