function CBIG_CreateLooseGCAHeadMask(input_file, output_file, dist)

% Create loose GCA head mask. This function 1) set nonzero value to 1; 2)
% dilate the mask; 3) fill in the holes
%
%   CBIG_CreateLooseGCAHeadMask(input_file, output_file, dist)
%   Input:
%       input_file  : input mri nifti file
%       output_file : output mri nifti file
%       dist        : dilate how many voxels
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(dist))
   dist = str2num(dist); 
end

% set non background to be 0.
disp('Create initial mask');
mask = MRIread(input_file);
mask.vol(mask.vol~=0) = 1;
MRIwrite(mask, output_file);

% dilate mask
disp(['Dilate mask by ' num2str(dist)]);
CBIG_DilateMask(output_file, output_file, dist)

% fill hole
disp('Fill holes');
mask = MRIread(output_file);
mask.vol = imfill(mask.vol);
MRIwrite(mask, output_file);
system(['mri_binarize --match 1 --i ' output_file ' --o ' output_file]);

exit
