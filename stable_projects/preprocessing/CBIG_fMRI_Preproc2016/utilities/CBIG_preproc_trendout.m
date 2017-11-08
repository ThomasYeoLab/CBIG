function output = CBIG_preproc_trendout(Y)
%This function fit a linear line for each time course (y=ax+b+e). Then, remove the mean
%and linear slope and get the resid (resid = y-ax-b).
%Y is a TxN matrix, T is num of time points, N is num of voxels. 
%
%  output = CBIG_preproc_trendout(Y)
%  Example:mri=MRIread('Sub0015.nii.gz'); 
%          vol=mri.vol; vol_size=size(vol)
%          Y=reshape(mri.vol,[vol_size(1)*vol_size(2)*vol_size(3) vol_size(4)])';
%          output = CBIG_preproc_trendout(Y)
%
%  Author: Nanbo Sun 
%  Date: 2016/06/13
%  
%  PS: The function of this code is a replication of trendout.awk from Washington University. 
%  The trendout.awk assumes that input Y is zero mean.


[n,ncol] = size(Y);

% get a -1 to 1 array with the same time points as input
x = zeros(n,1);
for i=0:(n-1)
    x(i+1) = -1+2*i/(n-1);
end

% fit a line and get sxx=sum((x_i)^2) sxy=sum(x_i*Y_i)
sxx = n*(n+1)/(3*(n-1));
sy = sum(Y);
sxy = sum(bsxfun(@times,Y,x));

% get the slope
a0 = sy/n;
a1 = sxy/sxx;

% regress out the line and get the residual
output = Y - x*a1;
