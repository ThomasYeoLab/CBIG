function CBIG_UpsampleVolume(input, output, scales, interp_type)

% Upsample a volume.
% 
%   CBIG_UpsampleVolume(input, output, scales, interp_type)
%   Input:
%       input       : input nifti file
%       output      : output nifti file
%       scales      : scales = FS_vox (order as seen on freeview)
%       interp_type : interpolation cubic, trilin, nearest (def is trilin)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(length(scales) ~= 3)
   error('Function assumes 3D scaling'); 
end

if(mod(scales(1), 2) == 0)
    warning('Requires odd scale for scaling to be centered');
end

if(mod(scales(2), 2) == 0)
    warning('Requires odd scale for scaling to be centered');
end

if(mod(scales(3), 2) == 0)
    warning('Requires odd scale for scaling to be centered');
end

x = MRIread(input);
sizes = size(x.vol);

x.vol = zeros(sizes .* [scales(2) scales(1) scales(3)]);
x.volres = x.volres ./ scales;

x.vox2ras(1:3, 1) = x.vox2ras(1:3, 1) / scales(1);
x.vox2ras(1:3, 2) = x.vox2ras(1:3, 2) / scales(2);
x.vox2ras(1:3, 3) = x.vox2ras(1:3, 3) / scales(3);

x.vox2ras0(1:3, 1) = x.vox2ras0(1:3, 1) / scales(1);
x.vox2ras0(1:3, 2) = x.vox2ras0(1:3, 2) / scales(2);
x.vox2ras0(1:3, 3) = x.vox2ras0(1:3, 3) / scales(3);

x.vox2ras1(1:3, 1) = x.vox2ras1(1:3, 1) / scales(1);
x.vox2ras1(1:3, 2) = x.vox2ras1(1:3, 2) / scales(2);
x.vox2ras1(1:3, 3) = x.vox2ras1(1:3, 3) / scales(3);

MRIwrite(x, output);

system(['mri_vol2vol --mov ' input ' --targ ' output ' --o ' output ' --regheader --interp ' interp_type]);