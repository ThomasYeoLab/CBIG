function CBIG_IndCBM_write_cerebellum_dlabel(surf_labels, cbm_labels, template_file, colortable, output_name)

% CBIG_IndCBM_write_cerebellum_dlabel(surf_labels, cbm_labels, template_file, colortable, output_name)
%
% This function write parcellation labels into a dlabel file. Including the
% cerebral cortical parcellation and the cerebellar parcellation. The
% parcellation will be colored based on the given color table. 
%
% Input:
%     - surf_labels: 
%           2N x 1 vector. Parcellation labels of the cereberal cortex.
%           Left hemisphere and right hemisphere are concatenated.
% 
%     - cbm_labels: 
%           M x 1 vector. Parcellation labels of the cerebellum.
%
%     - template_file: 
%           Cifti template dscalar file for a given mesh and specifying the
%           cerebellar voxels. Please create this file before you run this 
%           function. See CBIG_IndCBM_create_template.sh for how to create 
%           this template file. 
%           Example:
%           './examples/example_files/Sub1_fsaverage5_cerebellum_template.dscalar.nii'
%
%     - colortable: 
%           Struct with 'struct_names' and 'table'. 
%           'struct_names{i}' is the name of label i.
%           Each row of 'table' is the RGB values of label i.
%
%     - output_name: 
%           Specify the output file name. 
%
% Example:
% CBIG_IndCBM_write_cerebellum_dlabel([lh_labels; rh_labels], cbm_labels, 'proj/template.dscalar.nii', colortable, ...
% 'Sub1_parcellation')
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cifti = ft_read_cifti(template_file);
k = length(colortable.struct_names);

% create a file of label structure name and RGB color
fid = fopen([output_name,'_info.txt'], 'w');
for i = 1:k
    fprintf(fid, '%s\n', colortable.struct_names{i});
    fprintf(fid, '%i %i %i %i %i\n', i, colortable.table(i, 1), colortable.table(i, 2), colortable.table(i, 3), 255);
end
fclose(fid);

cifti.dscalar = double([surf_labels; cbm_labels]);

% write into a cifti file
ft_write_cifti(output_name, cifti, 'parameter', 'dscalar');

% apply wb_commnand to generate .dlabel file
input = [output_name, '.dscalar.nii'];
label_list_file = [output_name,'_info.txt'];
output = [output_name, '.dlabel.nii'];

command = ['wb_command -cifti-label-import ', input, ' ', label_list_file, ' ', output];
system(command);
delete(input, label_list_file);

end