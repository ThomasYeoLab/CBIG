function surf2vol_correlation_profile = CBIG_SmoothCorrelationProfileSurf2Vol2mm(surf2vol_correlation_profile, mask, method, filter_size, arg)

% surf2vol_correlation_profile = CBIG_SmoothCorrelationProfileSurf2Vol2mm(surf2vol_correlation_profile, mask, method, filter_size, arg)
% This function is used to smooth a surface to volume correlation profile
% with filter_size.
% surf2vol_correlation_profile = N x ROI
% method = SMOOTH => Everything outside the mask is smoothed ignoring the mask.
% method = SAME => Everything outside mask is set to same values as before.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


tmpVol = mask.vol;
mask_index = find(mask.vol == 1);
for i = 1:size(surf2vol_correlation_profile, 2)
   disp(num2str(i));
   tmpVol(mask_index) = surf2vol_correlation_profile(:, i);

   if(nargin == 5)
       tmpVol = CBIG_Smooth3DVolumeWithMask(tmpVol, mask.vol, method, filter_size, arg); 
   else
       tmpVol = CBIG_Smooth3DVolumeWithMask(tmpVol, mask.vol, method, filter_size);
   end
   
   surf2vol_correlation_profile(:, i) = tmpVol(mask_index);
end

