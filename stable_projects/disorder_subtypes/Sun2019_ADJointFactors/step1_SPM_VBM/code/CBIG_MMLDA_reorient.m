% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% reorient all images in the raw folder
spm_tool_dir = [spm_dir '/toolbox/OldNorm/']
addpath(spm_dir)
addpath(spm_tool_dir)

% roughly translate the image
x = load(reorient_sample);
f = strtrim(image_path);
N = nifti(f);
N.mat = x.M * N.mat;
create(N);

% reorient the image so that origin is in the AC accurately (Optional)
M_auto = CBIG_MMLDA_auto_reorient(image_path); 

% save the reorientation matrix 
M = M_auto*x.M;
[filepath, name, ext] = fileparts(image_path);
save([filepath '/' name '.mat'], 'M')

rmpath(spm_dir)
rmpath(spm_tool_dir)

