function CBIG_Smooth4DVolume(input_file, output_file, mask_file, outside_mask_type, method, filter_size, arg)

% CBIG_Smooth4DVolume(input_file, output_file, mask_file, outside_mask_type, method, filter_size, arg)
% This function is used to smooth a 4D data with filter_size.
% outside_mask_type = SMOOTH => Everything outside the mask is smoothed ignoring the mask.
% outside_mask_type = SAME => Everything outside mask is set to same values as before.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


input = MRIread(input_file);

if(strcmp(mask_file, 'NONE'))
    mask = [];
else
    mask = MRIread(mask_file);
end
    
if(ischar(filter_size))
   filter_size = str2num(filter_size); 
end

if(nargin >= 7)
    if(ischar(arg))
        arg = str2num(arg);
    end
end

disp('Smoothing...');
for i = 1:size(input.vol, 4)
    if(nargin < 7)
        if(isempty(mask))
            input.vol(:, :, :, i) = smooth3(squeeze(input.vol(:, :, :, i)), method, filter_size);
        else
            input.vol(:, :, :, i) = CBIG_Smooth3DVolumeWithMasks(squeeze(input.vol(:, :, :, i)), mask.vol, outside_mask_type, method, filter_size);
        end
    else
        if(isempty(mask))
            input.vol(:, :, :, i) = smooth3(squeeze(input.vol(:, :, :, i)), method, filter_size, arg);
        else
            input.vol(:, :, :, i) = CBIG_Smooth3DVolumeWithMasks(squeeze(input.vol(:, :, :, i)), mask.vol, outside_mask_type, method, filter_size, arg);
        end
    end
end

MRIwrite(input, output_file);

