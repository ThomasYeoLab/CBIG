% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% reorient all images in the raw folder
spm_tool_dir = [spm_dir '/toolbox/OldNorm/']
addpath(spm_dir)
addpath(spm_tool_dir)

% apply reorient matrix to T1 images
x = load(reorient_matrix);
f = strtrim(image_path);
N = nifti(f);
N.mat = x.M * N.mat;
create(N);

rmpath(spm_dir)
rmpath(spm_tool_dir)

