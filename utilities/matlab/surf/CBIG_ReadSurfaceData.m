function [surf_data, data_struct, data_size] = CBIG_ReadSurfaceData(input_file)

% [surf_data, data_struct, data_size] = CBIG_ReadSurfaceData(input_file)
%
% This function reads the content of surface data, including NIFTI and
% GIFTI format. For NIFTI file, the data in 'vol' can be either 3D or 4D.
% If it is 3D (x,y,z), it will be reshaped to (x*y*z,1). If it is 4D (x,y,z,t), 
% it will be reshaped to (x*y*z,t). No reshaping is performed for GIFTI
% file.
% 
% Input:
%     - input_file:
%       The full path of input file name. It can be in either NIFTI
%       (.nii) or GIFTI (.gii) format.
%
% Output:
%     - surf_data:
%       A num_voxels x num_timepoints matrix which is the content of
%       data_struct.vol (for NIFTI) or data_struct.cdata (for GIFTI).
%
%     - data_struct:
%       The structure read in by MRIread() or gifti().
%
%     - data_size:
%       The size of data_struct.vol (NIFTI) before reshaping or data_struct.cdata (GIFTI).
%
% Example:
%     - NIFTI input:
%       [surf_data, data_struct, data_size] = CBIG_ReadSurfaceData('/path/to/file/input.nii.gz')
%
%     - GIFTI input:
%       [surf_data, data_struct, data_size] = CBIG_ReadSurfaceData('/path/to/file/input.gii')
%
% Written by Yan Hongwei and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if (contains(input_file, '.nii'))
    % if input_file is NIFTI file
    data_struct = MRIread(input_file);
    surf_data = data_struct.vol;
    data_size = size(surf_data);
    if(length(data_size) == 3)
        % 3D data (x,y,z), reshape to (x*y*z, 1)
        surf_data = reshape(surf_data, prod(data_size(1:3)), 1);
    elseif (length(data_size) == 4)
        % 4D data (x,y,z,t), reshape to (x*y*z, t)
        surf_data = reshape(surf_data, prod(data_size(1:3)), data_size(4));
    else
        error('Only 3D and 4D NIFTI data structures are supported');
    end
elseif (contains(input_file, '.gii'))
    % if input_file is GIFTI file
    data_struct = gifti(input_file);
    if (isfield(data_struct, 'cdata') == 0)
        warning('The field ''cdata'' is not detected. You should use data_struct for further processing');
        surf_data = [];
        data_size = [];
        return;
    end
    surf_data = data_struct.cdata;
    data_size = size(surf_data);
else
    % if input_file format is unrecognized.
    error('Input argument ''input_file'' shoule be NIFTI or GIFTI file.');
end

end