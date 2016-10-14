function CBIG_AnalyzeDirectAndIndirectCorrespondence

% CBIG_AnalyzeDirectAndIndirectCorrespondence
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

indirect_count = MRIread('1000sub.FSL_MNI152.1mm.count_map.nii.gz');
indirect_vertex = MRIread('1000sub.FSL_MNI152.1mm.vertex_map.nii.gz');

direct_count = MRIread('1000sub.FSL_MNI152_FS.DIRECT.count_map.nii.gz');
direct_vertex = MRIread('1000sub.FSL_MNI152_FS.DIRECT.vertex_map.nii.gz');

index = find(indirect_count.vol(:) > 0 & direct_count.vol(:) > 0);

lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');

indirect_vertices = zeros(3, length(index));
tmp = indirect_vertex.vol(index);
indirect_vertices(:, tmp > 0) = lh_avg_mesh.vertices(:, tmp(tmp > 0));
indirect_vertices(:, tmp < 0) = rh_avg_mesh.vertices(:, -tmp(tmp < 0));

direct_vertices = zeros(3, length(index));
tmp = direct_vertex.vol(index);
direct_vertices(:, tmp > 0) = lh_avg_mesh.vertices(:, tmp(tmp > 0));
direct_vertices(:, tmp < 0) = rh_avg_mesh.vertices(:, -tmp(tmp < 0));

diff = sqrt(sum((indirect_vertices - direct_vertices).^2, 1));
output = direct_vertex;
output.vol(:) = 0;
output.vol(index) = diff;
MRIwrite(output, 'diff.nii.gz');
