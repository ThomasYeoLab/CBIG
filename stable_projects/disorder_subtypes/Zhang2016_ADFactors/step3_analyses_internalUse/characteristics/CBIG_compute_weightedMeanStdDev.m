function [wMean, wStdDev] = CBIG_compute_weightedMeanStdDev(X, y)

% [wMean, wStdDev] = CBIG_compute_weightedMeanStdDev(X, y)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

no_subtypes = size(X, 2);


%--------- Compute weighted mean
wMean = zeros(1, no_subtypes);
for idx = 1:no_subtypes
    w = X(:, idx);
    a = y.*w;
    wMean(idx) = sum(a)/sum(w);
end

%--------- Compute weighted standard deviationss
wStdDev = zeros(1, no_subtypes);
for idx = 1:no_subtypes
    w = X(:, idx);
    a = (y-wMean(idx)).^2;
    b = sum(a.*w);
    c = (nnz(w)-1)/nnz(w)*sum(w);
    wStdDev(idx) = sqrt(b/c);
end
