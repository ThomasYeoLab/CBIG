function C = CBIG_preproc_corr_matrix(A,B)

% C = CBIG_preproc_corr_matrix(A,B)
%
% Given two matrixs, calculate the Pearson correlation between each corresponding column.
% Input:
%   - A: a TxN matrix, T is num of time points, N is num of voxels  
%   - B: a TxN matrix, T is num of time points, N is num of voxels  
% Output:
%   - C: a 1xN matrix, N is num of voxels
%
% Author: Nanbo Sun
% Date: 2016/06/15

% Pearson correlation
% r = sum((x_i-mean(x))(y_i-mean(y)))/(sqrt(sum((x_i-mean(x))^2))sqrt(sum((y_i-mean(y))^2)))

% demean these two matrixs
A=bsxfun(@minus,A,mean(A,1));
B=bsxfun(@minus,B,mean(B,1));

% follow the equation above
A=bsxfun(@times,A,1./sqrt(sum(A.^2,1)));
B=bsxfun(@times,B,1./sqrt(sum(B.^2,1)));
C=sum(A.*B,1);