function newPhi = CBIG_AuthorTopic_NormalizePhi(phi, params)
% newPhi = CBIG_AuthorTopic_NormalizePhi(phi, params)
%
% Normalize the parameters phi of the variational distribution.
%
% Input:
%  - phi   : current estimate of the parameters of the variational
%            distribution.
%  - params: struct containing parameters of the author-topic model and
%            the Collapsed Variational Bayes (CVB) algorithm.
% Output:
%  - newPhi: updated parameters of the variational distribution of the
%            CVB algorithm.
%
% Example:
%   newPhi = CBIG_AuthorTopic_NormalizePhi(phi, params)
%   Return the normalized parameters of the variational distribution.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  newPhi = cell(params.E, 1);
  for e = 1:params.E
    minPhi  = min(min(phi{e},[], 3), [], 2);
    maxPhi  = max(max(phi{e},[], 3), [], 2);
    meanPhi = mean(mean(phi{e}, 3), 2);
    flag    = (minPhi > 0) | (maxPhi < 0);

    newPhi{e} = phi{e};
    newPhi{e}(flag, :, :) = bsxfun(@minus, phi{e}(flag, :, :), meanPhi(flag));

    newPhi{e} = exp(newPhi{e});
    denom = sum(sum(newPhi{e}, 3), 2);
    tmp = exp(phi{e});
    newPhi{e}(isinf(denom), :, :) = tmp(isinf(denom), :, :);
    denom = sum(sum(newPhi{e}, 3), 2);
    newPhi{e} = bsxfun(@times, newPhi{e}, 1./ denom);
    newPhi{e}(denom == 0, :, :) = 0;
    if sum(isnan(newPhi{e}(:))) > 0
      disp('Break');
    end;
  end;
end
