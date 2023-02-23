function cov_xy = CBIG_PFM_cov_matrix(x,y)

% cov_xy = CBIG_PFM_cov_matrix(x,y)
%
% This function calculate covariance between x and y
% 
% Inputs:
%   - x
%     N*K matrix, N is # subjects, K is # features
%
%   - y
%     N*1 vector, N is # subjects
% 
% Outputs:
%   - cov_xy
%     K*1 covariance matrix
%
% Written by Kong Ru and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
x_demean = bsxfun(@minus,x,CBIG_nanmean(x));
x_demean_nan = ~isnan(x_demean);
x_demean(isnan(x_demean)) = 0;
cov_xy = x_demean'*bsxfun(@minus,y,mean(y))./(sum(x_demean_nan,1)-1)';
end

