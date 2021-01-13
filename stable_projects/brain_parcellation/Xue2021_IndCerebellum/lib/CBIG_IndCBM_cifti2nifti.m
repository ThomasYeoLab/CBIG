function CBIG_IndCBM_cifti2nifti(cifti_file, mask_file, output_file, colortable)

% CBIG_IndCBM_cifti2nifti(cifti_file, mask_file, output_file, colortable)
%
% This function write cerebellar parcellation/confidence in cifti file into
% a nifti file. 
%
% Input:
%     - cifti_file: 
%           Path of the cifti file needed to be converted. 
%           Example: ../examples/ref_output/sub1/IndCBM_parcellation_top100.dlabel.nii
% 
%     - mask_file: 
%           Path of the cerebellar mask file.
%           Example: ../examples/input/mask/sub1/sub1_bin_mask_4mm.nii.gz
%
%     - output_file: 
%           Path of output file. Could be <name>.nii or <name>.nii.gz
%
% Optional input:
%
%     - colortable: 
%           Struct with 'struct_names' and 'table'. 
%           'struct_names{i}' is the name of label i.
%           Each row of 'table' is the RGB values of label i.
%           If 'colortable' is given, this function will save a txt look up
%           table file for freeview. 
%
% Example:
% CBIG_IndCBM_cifti2nifti('proj/parcellation.dlabel.nii', 'proj/mask.nii.gz', 'proj/parcellation.nii.gz')
% CBIG_IndCBM_cifti2nifti('proj/parcellation.dlabel.nii', 'proj/mask.nii.gz', 'proj/parcellation.nii.gz', colortable)
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Load input files
cifti = ft_read_cifti(cifti_file);
mask = MRIread(mask_file);
output = mask;
output.vol = zeros(size(output.vol));

% Find cerebellar structures
for i=1:length(cifti.brainstructurelabel)
    if(strcmp(cifti.brainstructurelabel{i}, 'CEREBELLUM_LEFT'))
        cbm(1) = i;
    elseif(strcmp(cifti.brainstructurelabel{i}, 'CEREBELLUM_RIGHT'))
        cbm(2) = i;
    end
end

% Assign labels in the volume
index = (cifti.brainstructure == cbm(1)) | (cifti.brainstructure == cbm(2));
label = cifti.dscalar(index);
ras_list = cifti.pos(index,:);
N = length(label);
for i = 1:N
    ras = ras_list(i,:)';
    vox = CBIG_ConvertRas2Vox(ras, cifti.transform);
    vox = ceil(vox([2 1 3]));
    % Note that cifti.pos already add one voxel to RAS so there's no need
    % to add 1 after converted to vox.
    output.vol(vox(1), vox(2), vox(3)) = label(i);
end

% Write to file
MRIwrite(output, output_file);

% Write look up table if colortable is given
if(nargin==4)
    filepath = fileparts(output_file);
    LUT_file = fullfile(filepath, 'Freesurfer_cerebellum_LUT.txt');
    if(exist(LUT_file, 'file'))
        delete(LUT_file)
    end
    fid = fopen(LUT_file, 'w');
    for i = 1:length(colortable.struct_names)
        fprintf(fid, '%g ', i);
        fprintf(fid, colortable.struct_names{i});
        fprintf(fid, ' %g %g %g 0\n', colortable.table(i,1), colortable.table(i,2), colortable.table(i,3));
    end
    fclose(fid);
end

end