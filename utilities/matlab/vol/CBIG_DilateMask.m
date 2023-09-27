function CBIG_DilateMask(input_file, output_file, dist, exit_flag)

% Dilate a binary mask.
%
%     CBIG_DilateMask(input_file, output_file, dist, exit_flag)
%     Input:
%         input_file    : input binary mask
%         output_file   : output file name
%         dist          : threshold for Euclidean distance after bwdist
%         exit_flag     : 1-exit; 0-not exit
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(nargin < 4)
   exit_flag = 0; 
end

if(ischar(dist))
   dist = str2num(dist); 
end

mask = MRIread(input_file);
dt = bwdist(mask.vol);
mask.vol(dt <= dist) = 1;
MRIwrite(mask, output_file);
system(['mri_binarize --match 1 --i ' output_file ' --o ' output_file]);

if(exit_flag)
    exit
end
