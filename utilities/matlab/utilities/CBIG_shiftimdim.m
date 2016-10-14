function CBIG_shiftimdim(input_file, output_file, shiftdim1, shiftdim2, shiftdim3)

% Circle shift a 3D nifti file.
% 
%   CBIG_shiftimdim(input_file, output_file, shiftdim1, shiftdim2, shiftdim3)
%   Input:
%       input_file  : input file
%       output_file : output file
%       shiftdim1   : shift number in dimension 1
%       shiftdim2   : shift number in dimension 2
%       shiftdim2   : shift number in dimension 3
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


shiftdim1 = str2num(shiftdim1);
shiftdim2 = str2num(shiftdim2);
shiftdim3 = str2num(shiftdim3);

input = MRIread(input_file);
input.vol = circshift(input.vol, [shiftdim1 shiftdim2 shiftdim3]);
MRIwrite(input, output_file);

exit
