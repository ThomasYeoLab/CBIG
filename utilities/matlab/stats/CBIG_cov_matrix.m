function cov_xy = CBIG_cov_matrix(x,y)

% function cov_xy = CBIG_cov_matrix(x,y)
% 
% Calculates the covariance matrix given two matrices x and y.
%
% Inputs:
% - x 
%   N*K1 matrix, N is # observation, K1 is # variables
%
% - y
%   N*K2 matrix, N is # observation, K2 is # variables
%
% Outputs:
% - cov_xy: 
%   K1*K2 covariance matrix
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cov_xy = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);

end