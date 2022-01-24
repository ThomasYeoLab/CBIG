function cov_xy = CBIG_TRBPC_compute_covariance(x,y)

% cov_xy = CBIG_TRBPC_compute_covariance(x,y)
%
% x: N*K1 matrix, N is # observation, K is # variables
% y: N*K2 matrix, N is # observation, K is # variables
% cov_xy: K1*K2 covariance matrix
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cov_xy = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);

end