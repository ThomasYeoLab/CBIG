function [] = CBIG_write_correlation_matrix_into_ciftiformat(corrmat, output_name, output_dir)
%function [] = CBIG_write_correlation_matrix_into_ciftiformat(corrmat, output_name, output_dir)
%
% This function will take in a correlation matrix (.mat file) and convert
% it into cifti format (.dconn.nii file). 
%
% After getting the .dconn.nii file, you may load it into 'wb_view' for 
% interactive seed selection and analysis.
% 
% To visulize the correlation matrix in wb_view, you also need two surface
% files as underlay:
% 1) fsLR space
%    Surfaces files can be found under '$CBIG_CODE_DIR/data/templates/surface'
%    You can choose the surface files according to your need, i.e.
%    'fs_LR_32k/fsaverage.<?>.very_inflated.32k_fs_LR.surf.gii'
%
% 2) fsaverage space
%   You need to create two '.surf.gii' files that 'wb_view' can take in.
%   You can find the fsaverage surf files you need under this path:
%   '$FREESURFER_HOME/subjects/subjects/fsaverage<?>/surf/'
%   Then you can convert them into 'surf.gii' format using 'mris_convert' command 
%   and load them into 'wb_view' as underlay.
%
%
% Input:
%       corrmat     = your input NxN correlation matrix
%                     a) matrix should have both lh and rh (i.e. [lh2lh, lh2rh; rh2lh, rh2rh])
%                     b) matrix can be in fsaverage space or fsLR space
%                     c) medial wall should not be removed in the matrix
%       output_name = name of the output file
%       output_dir  = output directory
%
% Output:
%       a cifti file: '<output_name>.dconn.nii' saved under output_dir
%
% Written by Yang-Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%% construct a cifti template
resolution = size(corrmat, 1) / 2;

% load in a gifti template
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
g = gifti(fullfile(CBIG_CODE_DIR, 'data', 'templates', 'surface', 'fs_LR_32k', ...
    'example_func', 'Parcels_L.func.gii'));

% construct two gifti metric files
metric = single(ones(resolution, 1));
g.cdata = metric;
save(g, fullfile(output_dir, 'L.example.func.gii'));
save(g, fullfile(output_dir, 'R.example.func.gii'));

% construct cifti template using wb_command
cifti_out = fullfile(output_dir, 'cifti_template.dscalar.nii');
left_metric = fullfile(output_dir, 'L.example.func.gii');
right_metric = fullfile(output_dir, 'R.example.func.gii');

command = ['wb_command -cifti-create-dense-scalar ' cifti_out ...
    ' -left-metric ' left_metric ' -right-metric ' right_metric];
system(command);


%% rewrite cifti template
% open cifti template
tmpfile = tempname;
command = ['wb_command -cifti-convert -to-gifti-ext ' cifti_out ' ' tmpfile '.gii'];
system(command);
cifti = gifti([tmpfile '.gii']);
% delete intermedaite files
delete([tmpfile '.gii'],[tmpfile '.gii.data']);

% replace data 
corrmat = single(corrmat);
cifti.cdata = corrmat;

% save into .dtseries.nii format
dtseries_file = [output_dir, '/', output_name, '.dtseries.nii'];
save(cifti, [dtseries_file '.gii'], 'ExternalFileBinary')
command = ['wb_command -cifti-convert -from-gifti-ext ' dtseries_file '.gii ' ...
    dtseries_file ' -reset-timepoints 1 0'];
system(command);
% delete intermedaite files
delete([dtseries_file '.gii'],[dtseries_file '.dat']);

clear corrmat


%% convert into .dconn.nii format
dconn_file = [output_dir, '/', output_name, '.dconn.nii'];

command = ['wb_command -cifti-copy-mapping ' dtseries_file ' ROW '...
    dtseries_file ' COLUMN ' dconn_file];
system(command);

% remove intermediate files
delete(left_metric, right_metric, cifti_out, dtseries_file);

end