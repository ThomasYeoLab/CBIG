function CBIG_RF_make_xyzIndex_vol(template, output_dir, output_prefix)

% CBIG_RF_make_xyzIndex_vol(mask, output_dir, output_prefix)
%
% This function creates index volumes from a template (e.g. MNI152)
% The index files contains matrix of x/y/z RAS coordinates of each voxel 
% in the template's space, i.e. with the same dimensions as the template
%
% Input:
%     - template     :
%                      absolute/relative path to the template file, 
%                      which should be a volumetric file that can be read by MRIread()
%     - output_dir   :
%                      absolute/relative path to directory where output should be stored
%     - output_prefix:
%                      desired prefix for the outputs
%
% Output:
%     - There is no function output.
%     - 3 index files are created in output_dir:
%           [output_prefix]_x.INDEX.nii.gz
%           [output_prefix]_y.INDEX.nii.gz
%           [output_prefix]_z.INDEX.nii.gz
%
% Example:
% CBIG_RF_make_xyzIndex_vol('~/templates/MNI152_T1_1mm_brain.nii.gz', '../index_MNI/', 'MNI_1mm')
% This command generates index volumes based on header information of MNI152_T1_1mm_brain.nii.gz. 
% The index volumes are placed in ../results/index_MNI folder and have prefix of 'MNI_1mm'
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 3
    disp('usage: CBIG_RF_make_xyzIndex_vol(mask, output_dir, output_prefix)');
    return
end

%Load volume and extract info
template = MRIread(template);
dimensions = size(template.vol);
indices_total = numel(template.vol);
vox2ras = template.vox2ras;

%Create the index matrix in Matlab matrix format
[mat_j, mat_i, mat_k] = meshgrid(1:dimensions(2), 1:dimensions(1), 1:dimensions(3));

%Convert the Matlab format to voxel format
vox_x = mat_j - 1;
vox_y = mat_i - 1;
vox_z = mat_k - 1;

%Convert the voxel format to RAS format
vox_x = reshape(vox_x, 1, indices_total);
vox_y = reshape(vox_y, 1, indices_total);
vox_z = reshape(vox_z, 1, indices_total);
ras = CBIG_ConvertVox2Ras([vox_x; vox_y; vox_z], vox2ras);
ras_x = reshape(ras(1, :), dimensions(1), dimensions(2), dimensions(3));
ras_y = reshape(ras(2, :), dimensions(1), dimensions(2), dimensions(3));
ras_z = reshape(ras(3, :), dimensions(1), dimensions(2), dimensions(3));

%Save the 3 index volumes separately
template.vol = ras_x;
MRIwrite(template, [output_dir '/' output_prefix '_x.INDEX.nii.gz']);
template.vol = ras_y;
MRIwrite(template, [output_dir '/' output_prefix '_y.INDEX.nii.gz']);
template.vol = ras_z;
MRIwrite(template, [output_dir '/' output_prefix '_z.INDEX.nii.gz']);

end
