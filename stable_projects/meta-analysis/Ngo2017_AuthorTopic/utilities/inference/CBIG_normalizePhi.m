function newPhi = normalizePhi(phi, params)
  % newPhi = normalizePhi(phi, params)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  newPhi = cell(params.D, 1);
  for d = 1:params.D
    minPhi  = min(min(phi{d},[],3),[],2);
    maxPhi  = max(max(phi{d},[],3),[],2);
    meanPhi = mean(mean(phi{d}, 3), 2);
    flag    = (minPhi > 0) | (maxPhi < 0);

    newPhi{d} = phi{d};
    newPhi{d}(flag,:,:) = bsxfun(@minus, phi{d}(flag,:,:), meanPhi(flag));
      
    newPhi{d} = exp(newPhi{d});
    denom = sum(sum(newPhi{d},3),2);
    tmp = exp(phi{d});
    newPhi{d}(isinf(denom),:,:) = tmp(isinf(denom),:,:);
    denom = sum(sum(newPhi{d},3),2);
    newPhi{d} = bsxfun(@times, newPhi{d}, 1./ denom);
    newPhi{d}(denom == 0,:,:) = 0;
    if sum(isnan(newPhi{d}(:))) > 0
      disp('Break');
    end;
  end;
end
