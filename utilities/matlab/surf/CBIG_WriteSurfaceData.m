function CBIG_WriteSurfaceData(output_file, surf_data, data_struct, data_size)

% CBIG_WriteSurfaceData(output_file, surf_data, data_struct, data_size)
%
% This function writes surface data to a file.
% 
% Input:
%     - output_file:
%       The output file name (full path). The output file should be either
%       in NIFTI (.nii) format or in GIFTI (.gii) format.
%
%     - surf_data:
%       The content of output signals. This should be in the shape of (n,t)
%       such that each row is surface data for one frame.
%
%     - data_struct:
%       The structure used for MRIwrite() if the output file is in NIFTI
%       format. 
%       The structure is ignored if the output file is in GIFTI format.
%
%     - data_size:
%       The size that the surf_data will be reshaped to for NIFTI format.
%       If it is not provided, the function will attempt to reshape
%       the data using data_struct.volsize.
%       This is ignored if the output file is in GIFTI format.
%
% Example:
%     - NIFTI file:
%       CBIG_WriteSurfaceData('/path/to/file/output.nii.gz', surf_data, data_struct, data_size)
%
%     - GIFTI file:
%       CBIG_WriteSurfaceData('/path/to/file/output.gii', surf_data)
%
% Written by Yan Hongwei and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(contains(output_file, '.nii'))
    % if output_file is NIFTI file    
    if (nargin < 3)
        error('Output format is NIFTI, ''data_struct'' is required argument');
    end
    if (isfield(data_struct, 'volsize') == 0)
        error('Output format is NIFTI, but ''data_struct'' is not NIFTI structure');
    end
    struct_data_size = data_struct.volsize;
    input_data_size = size(surf_data);
    if (nargin == 3)
        % if data_size not given, use data_struct.volsize to reshape        
        if (prod(struct_data_size) ~= input_data_size(1))
            error('data_size not given but data_struct.volsize does not match size of surf_data');
        end
        data_size = [struct_data_size, input_data_size(2)];
    else
        % if data_size is given, check consistency among data_size,
        % data_struct.volsize and surf_data
        len = length(data_size);
        if (prod(data_size) == prod(input_data_size))
            % if data_size is consistent with surf_data, use data_size for
            % reshape, check with data_struct.volsize for possible warning
            if (prod(struct_data_size) ~= prod(data_size(1:len-1)))
                warning('data_size and data_struct.volsize are inconsistent. Using data_size to reshape');
            end        
        else
            % if data_size is not consistent with surf_data, try to use
            % data_struct.volsize to reshape
            if (prod(struct_data_size) ~= input_data_size(1))
                error('both data_size and data_struct.volsize do not match size of surf_data');
            end
            warning('data_size and surf_data do not match. Using data_struct.volsize to reshape')
            data_size = [struct_data_size, input_data_size(2)];
        end
    end
    surf_data = reshape(surf_data, data_size);
    data_struct.vol = surf_data;
    MRIwrite(data_struct, output_file);
elseif (contains(output_file, '.gii'))
    % if output_file is GIFTI file
    data_struct = gifti(surf_data);
    save(data_struct, output_file);
else
    % if output_file format is unsupported
    error('Input argument ''output_file'' shoule be NIFTI or GIFTI file.');
end

end