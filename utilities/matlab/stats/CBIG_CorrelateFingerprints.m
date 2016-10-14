function output = CBIG_CorrelateFingerprints(mat1, mat2)

% output = CBIG_CorrelateFingerprints(mat1, mat2)
% Compute Pearson correlation between two matrixs column by column.
% 
%   CBIG_CorrelateFingerprints(mat1, mat2)
%   Input:
%       mat1    : Targets x seeds matrix
%       mat2    : Targets x seeds matrix
%   Output:
%       output  : 1 x seeds
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



mat1 = bsxfun(@minus, mat1, mean(mat1, 1));
mat1 = bsxfun(@times, mat1, 1./sqrt(sum(mat1.^2, 1)));

mat2 = bsxfun(@minus, mat2, mean(mat2, 1));
mat2 = bsxfun(@times, mat2, 1./sqrt(sum(mat2.^2, 1)));

output = sum(mat1 .* mat2, 1);