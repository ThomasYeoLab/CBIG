function CBIG_RF_example_compare(file1, file2)
% This function compares two nifti files and returns 1 if their .vol are completely the same
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

mri1 = MRIread(file1);
mri2 = MRIread(file2);
diff = sum(abs(mri1.vol(:) - mri2.vol(:)));

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
if diff > 0
  disp('%%  The two volumes are different!  %%')
else
  disp('%%  The two volumes are identical.  %%')
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
