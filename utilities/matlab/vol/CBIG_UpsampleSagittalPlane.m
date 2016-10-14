function CBIG_UpsampleSagittalPlane(input, output, scale, interp_type)

% Upsample a volume along sagittal plane.
% 
%   CBIG_UpsampleSagittalPlane(input, output, scale, interp_type)
%   Input:
%       input       : input nifti file
%       output      : output nifti file
%       scale       : a single parameter
%       interp_type : interpolation cubic, trilin, nearest (def is trilin)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(numel(scale) ~= 1)
   error('Scale has to be a single parameter'); 
end

if(mod(scale, 2) == 0)
    error('Requires odd scale for scaling to be centered');
end

x = MRIread(input);

if(sum(x.vox2ras(:) == 0) ~= 9)
    error('It is not possible to just upsample the sagittal plane because the voxels are not parallel to the ras axes');
end

% upsampling sagittal plane equivalent to upsampling anterior-posterior and
% superior-inferior axes

scale_vec = transpose(abs(scale * inv(x.vox2ras(1:3, 1:3)) * [0; 1; 1]));
scale_vec(scale_vec == 0) = 1;
CBIG_UpsampleVolume(input, output, scale_vec, interp_type); 