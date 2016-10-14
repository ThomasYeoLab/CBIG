function CBIG_GrabFrames(input_file, output_file, subsample)

% Extract several frames from a 4D nifti volume.
% 
%   CBIG_GrabFrames(input_file, output_file, subsample)
%   Input:
%       input_file  : input file name
%       output_file : output file name
%       subsample   : subsample vector (e.g. [1:10]) 
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(ischar(subsample))
   subsample = str2num(subsample); 
end

disp(['Extracting frames: ' num2str(subsample)]);

input = MRIread(input_file);
input.vol = input.vol(:, :, :, subsample);
input.nframes = length(subsample);

MRIwrite(input, output_file);


