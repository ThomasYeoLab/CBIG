function CBIG_ExtendValuesInMaskToEntireVolume(input_file, mask_file, output_file, exit_flag)

% Extend the values within the mask to the entire volume.
%
%     CBIG_ExtendValuesInMaskToEntireVolume(input_file, mask_file, output_file, exit_flag)
%     Input:
%         input_file    : input nifti file
%         mask_file     : mask nifti file
%         output_file   : output nifti file
%         exit_flag     : 1-exit matlab after the function; 0-doesn't exit matlab after the function 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(exit_flag))
    exit_flag = str2num(exit_flag); 
end

input = MRIread(input_file);
mask = MRIread(mask_file);

if(max(max(max(abs(size(mask.vol) - size(input.vol))))) ~= 0)
    error('input file and mask not the same dimensions'); 
end

[dist, id] = bwdist(mask.vol);
output = input;
output.vol = input.vol(id);
MRIwrite(output, output_file);

if(exit_flag)
    exit 
end

