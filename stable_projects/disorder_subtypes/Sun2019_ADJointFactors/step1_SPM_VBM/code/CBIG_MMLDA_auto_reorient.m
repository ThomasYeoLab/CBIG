function M = CBIG_MMLDA_auto_reorient(path) 
% This code is from https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;d1f675f1.0810

% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

spmDir = which('spm');
spmDir = spmDir(1:end-5);
tmpl = [spmDir 'canonical/avg152T1.nii'];
vg = spm_vol(tmpl);
flags.regtype = 'rigid';

f = strtrim(path);
[filepath, name, ext] = fileparts(f);
tempfile = [filepath '/' name '_temp.nii'];
spm_smooth(f,tempfile,[12 12 12]);
vf = spm_vol(tempfile);

% rigid tranformation from native space to template space
[M,scal] = spm_affreg(vg,vf,flags);
M3 = M(1:3,1:3);
[u s v] = svd(M3);
M3 = u*v';
M(1:3,1:3) = M3;
N = nifti(f);
N.mat = M*N.mat;
create(N);
