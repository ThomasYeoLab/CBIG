function CBIG_WriteToReferenceVolume(input_file, ref_file, output_file, exit_flag)

% Write input volume to reference nifti file.
% 
%   CBIG_WriteToReferenceVolume(input_file, ref_file, output_file, exit_flag)
%   Input:
%       input_file  : input file
%       ref_file    : reference file
%       output_file : output file
%       exit_flag   : 1-exit matlab; 0-not
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(nargin < 4)
   exit_flag = 0; 
end

input = MRIread(input_file);
ref = MRIread(ref_file);

if(size(ref.vol) ~= size(input.vol))
   error('Dimensions not the same!'); 
end
ref.vol = input.vol;
MRIwrite(ref, output_file);

if(exit_flag)
   exit; 
end
