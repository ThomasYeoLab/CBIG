function x = CBIG_StableAtanh(x)

% Stable version of atanh. It constrains the input within [-1 1] and output within [atanh(-1+eps) atanh(1-eps)].
%
%     x = CBIG_StableAtanh(x)
%     Input:
%         x: a matrix
%     Output:
%         x: a matrix
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


x(x > 1) = 1;
x(x < -1) = -1;
x = atanh(x);
x(isinf(x) & x > 0) = atanh(1-eps);
x(~isreal(x) & real(x) > 0) = atanh(1-eps);

x(isinf(x) & x < 0) = atanh(-1+eps);
x(~isreal(x) & real(x) < 0) = atanh(-1+eps);
