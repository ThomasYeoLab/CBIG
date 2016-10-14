function [sbb, sww, icc] = CBIG_manova_cov(data)

% [sbb, sww, icc] = CBIG_manova_cov(data)
% Implemented as in Equation 3 of Amemiya 1985 (What should be done when an estimated
% between-group covariance Matrix is not Nonnegative Definite)
% WITHOUT correction for negative-definite covariance
%
% assume data = n x r x p 
% where n is # targets, r is # judges, p is # dimensions of observations
%
% Model assumption
%
% Y_{ij} = mu + b_i + w_{ij},  1 <= j <= r, 1 <= i <= n
%
% b_i ~ mean 0, cov s_bb;  w_{ij} ~ mean 0, cov s_ww
%
% p can be equal to 1.
%
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



n = size(data, 1);
r = size(data, 2);

if(n < 2)
   error('n has to be at least 2'); 
end

if(r < 2)
   error('r has to be at least 2'); 
end

meanYi = mean(data, 2); % n x p
meanY = squeeze(mean(meanYi, 1)); % 1 x p


% mbb = 0;
% diffY = bsxfun(@minus, meanYi, meanY);
% for i = 1:n
%    mbb = mbb + diffY(i, :)'*diffY(i, :); 
% end
% mbb = mbb/(n - 1)*r ;
diffY = bsxfun(@minus, squeeze(meanYi), meanY);
mbb = diffY'*diffY/(n-1)*r;


diffY = bsxfun(@minus, data, meanYi);
mww = 0;
for i = 1:n
    for j = 1:r
        mww = mww + squeeze(diffY(i, j, :))'*squeeze(diffY(i, j, :));
    end
end
mww = mww/n/(r-1);

% outputs
sww = mww;
sbb = 1/r*(mbb - mww);

% icc
icc = det(sbb)/det(sww+sbb);

