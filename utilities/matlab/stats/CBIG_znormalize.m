function CBIG_znormalize(input_file, output_file, mask_file)

% Znormalize the input nifti file along 1st dimension.
% 
%   CBIG_znormalize(input_file, output_file, mask_file)
%   Input:
%       input_file  : input file
%       output_file : output file
%       mask_file   : mask file
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


x = MRIread(input_file);

if(nargin < 3)
    mask_index = 1:numel(x.vol);
else
    mask = MRIread(mask_file);
    mask_index = find(mask.vol~=0);
end

val = x.vol(mask_index);
val = (val - mean(val))./std(val);
x.vol(mask_index) = val;

if(nargin == 3)
    x.vol(mask.vol==0) = min(val); 
end

MRIwrite(x, output_file);

exit
