function smoothedV = CBIG_Smooth3DVolumeWithMasks(V, mask, outside_mask_type, method, filter_size, arg)

% smoothedV = CBIG_Smooth3DVolumeWithMasks(V, mask, outside_mask_type, method, filter_size, arg)
% outside_mask_type = SMOOTH => Everything outside the mask is smoothed ignoring the mask.
% outside_mask_type = SAME => Everything outside mask is set to same values as before.
% Everything outside the mask is smoothed ignoring the mask
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



vals = unique(mask);
if(1+max(vals) ~= length(vals))
    error('Assumes mask contains values from 0 to N');
end

smoothedV = zeros(size(V));
for i = 1:max(vals)
   
    tmpV = V; tmpV(mask ~= i) = 0; % set up volume with 0 outside of mask == i
    tmpMask = ones(size(mask)); tmpMask(mask ~= i) = 0; % set up binary mask with 0 outside of mask == i
    
    if(nargin < 6)
        smoothed_mask = smooth3(tmpMask, method, filter_size);
        smoothedTmpV = smooth3(tmpV, method, filter_size);
    else
        smoothed_mask = smooth3(tmpMask, method, filter_size, arg);
        smoothedTmpV = smooth3(tmpV, method, filter_size, arg);
    end
    
    % adjusts for lower values due to boundaries
    smoothedV(mask == i) = smoothedTmpV(mask == i)./smoothed_mask(mask == i);
end

if(strcmp(outside_mask_type, 'SAME'))
  smoothedV(mask == 0) = V(mask == 0);
elseif(strcmp(outside_mask_type, 'SMOOTH'))
  if(nargin < 6)
    outsideV = smooth3(V, method, filter_size);
  else
    outsideV = smooth3(V, method, filter_size, arg);
  end
  smoothedV(mask == 0) = outsideV(mask == 0);
else
  error(['Does not handle outside mask type: ' outside_mask_type]);
end




