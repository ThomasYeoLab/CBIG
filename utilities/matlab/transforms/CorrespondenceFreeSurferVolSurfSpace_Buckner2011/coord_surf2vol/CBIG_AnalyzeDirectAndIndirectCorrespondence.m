function CBIG_AnalyzeDirectAndIndirectCorrespondence

% CBIG_AnalyzeDirectAndIndirectCorrespondence
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

load lh.1000sub.FSL_MNI.ras.mat
indirect_ras = ras;
load lh.FSL_MNI.DIRECT.ras.mat
direct_ras = ras;
tmp = sqrt(sum((direct_ras - indirect_ras).^2, 1));
lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage', 'white', 'cortex');
TrisurfMeshData(lh_avg_mesh, tmp, 1); shading interp; colorbar



load rh.1000sub.FSL_MNI.ras.mat
indirect_ras = ras;
load rh.FSL_MNI.DIRECT.ras.mat
direct_ras = ras;
tmp = sqrt(sum((direct_ras - indirect_ras).^2, 1));
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage', 'white', 'cortex');
TrisurfMeshData(rh_avg_mesh, tmp, 1); shading interp; colorbar
