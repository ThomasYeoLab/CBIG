function CBIG_SmoothCorrelationProfileSurf2Vol2mmWrapper(surf2vol_correlation_profile_file, output_profile, mask_file, method, filter_size, arg)

% CBIG_SmoothCorrelationProfileSurf2Vol2mmWrapper(surf2vol_correlation_profile_file, output_profile, mask_file, method, filter_size, arg)
% This function is used to smooth a surface to volume correlation profile
% with filter_size.
% surf2vol_correlation_profile is a .mat file contains a N x ROI matrix
% method = SMOOTH => Everything outside the mask is smoothed ignoring the mask.
% method = SAME => Everything outside mask is set to same values as before.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ischar(filter_size))
   filter_size = str2num(filter_size); 
end

if(nargin == 6 && ischar(arg))
   arg = str2num(arg); 
end

mask = MRIread(mask_file);
load(surf2vol_correlation_profile_file);


smoothed_profile = zeros(size(surf2vol_correlation_profile));
tmpVol = mask.vol;
mask_index = find(mask.vol == 1);
for i = 1:size(surf2vol_correlation_profile, 2)
   disp(num2str(i));
   tmpVol(mask_index) = surf2vol_correlation_profile(:, i);

   if(nargin == 6)
       tmpVol = CBIG_Smooth3DVolumeWithMask(tmpVol, mask.vol, method, filter_size, arg); 
   else
       tmpVol = CBIG_Smooth3DVolumeWithMask(tmpVol, mask.vol, method, filter_size);
   end
   
   smoothed_profile(:, i) = tmpVol(mask_index);
end

if(nargin == 6)
    surf2vol_correlation_profile = CBIG_SmoothCorrelationProfileSurf2Vol2mm(surf2vol_correlation_profile, mask, method, filter_size, arg);
else
    surf2vol_correlation_profile = CBIG_SmoothCorrelationProfileSurf2Vol2mm(surf2vol_correlation_profile, mask, method, filter_size);
end
save(output_profile, 'surf2vol_correlation_profile');
exit
