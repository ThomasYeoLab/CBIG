function CBIG_Smooth3DFrames(input_prefix, start, stop, input_file_type, output_prefix, mask_file, outside_mask_type, method, filter_size, arg)

% CBIG_Smooth3DFrames(input_prefix, start, stop, input_file_type, output_prefix, mask_file, outside_mask_type, method, filter_size, arg)
% This function is used to smooth a 3D data with filter_size for frames
% between start and stop.
% outside_mask_type = SMOOTH => Everything outside the mask is smoothed ignoring the mask.
% outside_mask_type = SAME => Everything outside mask is set to same values as before.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(strcmp(mask_file, 'NONE'))
    mask = [];
else
    mask = MRIread(mask_file);
end
    
if(ischar(filter_size))
   filter_size = str2num(filter_size); 
end

if(nargin >= 10)
    if(ischar(arg))
        arg = str2num(arg);
    end
end

if(ischar(start))
   start = str2num(start); 
end

if(ischar(stop))
   stop = str2num(stop); 
end

disp('Smoothing indivdual frames: ');
for i = start:stop
    input_file = [input_prefix num2str(i, '%04d') input_file_type];
    input = MRIread(input_file);
    
    if(i ~= stop)
        fprintf([num2str(i) ',']);
    else
        fprintf([num2str(i) '\n']);
    end
    
    if(nargin < 10)
        if(isempty(mask))
            input.vol = smooth3(input.vol, method, filter_size);
        else
            input.vol = CBIG_Smooth3DVolumeWithMasks(input.vol, mask.vol, outside_mask_type, method, filter_size);
        end
    else
        if(isempty(mask))
            input.vol = smooth3(input.vol, method, filter_size, arg);
        else
            input.vol = CBIG_Smooth3DVolumeWithMasks(input.vol, mask.vol, outside_mask_type, method, filter_size, arg);
        end
    end
    
    output_name = [output_prefix num2str(i, '%04d') input_file_type];
    MRIwrite(input, output_name);
end
fprintf('Done!! \n');

