function CBIG_SPGrad_create_cifti_template(resolution, output_dir)

% CBIG_SPGrad_create_cifti_template(resolution, output_dir)
%
% This script creates a cifti template for a given resolution
%
% Input:
%     - resolution: (scalar)
%       the number of vertices of surface mesh, e.g. 40962.
%
%     - output_dir:
%       the output directory which saves the cifti template.
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% load in a gifti template
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
g = gifti([CBIG_CODE_DIR, '/data/templates/surface/fs_LR_32k/example_func/Parcels_L.func.gii']);

% construct two gifti metric files
metric = single(ones(resolution, 1));
g.cdata = metric;
save(g, [output_dir, '/L.example.func.gii']);
save(g, [output_dir, '/R.example.func.gii']);

% construct cifti template using wb_command
cifti_out = [output_dir '/cifti_template.dscalar.nii'];
left_metric = [output_dir, '/L.example.func.gii'];
right_metric = [output_dir, '/R.example.func.gii'];

command = ['wb_command -cifti-create-dense-scalar ' cifti_out ...
    ' -left-metric ' left_metric ' -right-metric ' right_metric];
system(command);