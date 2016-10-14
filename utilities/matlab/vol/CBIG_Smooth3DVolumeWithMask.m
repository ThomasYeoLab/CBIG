function smoothedV = CBIG_Smooth3DVolumeWithMask(V, mask, outside_mask_type, method, filter_size, arg)

% smoothedV = CBIG_Smooth3DVolumeWithMask(V, mask, outside_mask_type, method, filter_size, arg)
% outside_mask_type = SMOOTH => Everything outside the mask is smoothed ignoring the mask.
% outside_mask_type = SAME => Everything outside mask is set to same values as before.
% mask assumed to be binary.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



tmpV = V;
tmpV(mask ~= 1) = 0;
if(nargin < 6)
    smoothed_mask = smooth3(mask, method, filter_size);
    smoothedV = smooth3(tmpV, method, filter_size);
else
    smoothed_mask = smooth3(mask, method, filter_size, arg);
    smoothedV = smooth3(tmpV, method, filter_size, arg);
end


smoothedV(mask == 1) = smoothedV(mask == 1)./smoothed_mask(mask==1);
if(strcmp(outside_mask_type, 'SAME'))
  smoothedV(mask ~= 1) = V(mask ~= 1);
elseif(strcmp(outside_mask_type, 'SMOOTH'))
  if(nargin < 6)
    outsideV = smooth3(V, method, filter_size);
  else
    outsideV = smooth3(V, method, filter_size, arg);
  end
  smoothedV(mask ~= 1) = outsideV(mask~=1);
else
  error(['Does not handle outside mask type: ' outside_mask_type]);
end
