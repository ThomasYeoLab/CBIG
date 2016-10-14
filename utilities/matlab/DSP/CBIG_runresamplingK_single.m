function out = CBIG_runresamplingK_single(x,Kmax,Kmin,B,d,nruns,epsilon)

% out = CBIG_runresamplingK_single(x,Kmax,Kmin,B,d,nruns,epsilon)
%
% Random permute data and split vertices into halves. Perform clustering
% algorithm on each half.
%
% Input arguments:
%   - x       : input profiles data.
%   - Kmax    : maximum number of clusters.
%   - Kmin    : minimum number of clusters.
%   - B       : number of different random permutations.
%   - d       : the dimensionality. If enter zero, the actual dimensionality
%               of the data.
%   - nruns   : number of repetitions when using random initializations.
%   - epsilon : value of the condition for breaking the the iterations.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


rand('seed',100*sum(clock));


N = size(x,1);
Nhalf = ceil(N/2);

for k = Kmin:Kmax
    for i=1:B
        
        indrand = randperm(size(x,1));        
        
        res = direcClus_fix_bessel_bsxfun(x(indrand(1:Nhalf),:),k,d,nruns,0,0,0,epsilon,1);
        out{k,i,1}.inds = indrand(1:Nhalf);
        out{k,i,1}.p = res.p;
        out{k,i,1}.mtc = res.mtc;
        out{k,i,1}.lambda = res.lambda;
        
        res = direcClus_fix_bessel_bsxfun(x(indrand(Nhalf+1:end),:),k,d,nruns,0,0,0,epsilon,1);
        out{k,i,2}.inds = indrand(Nhalf+1:end);
        
        out{k,i,2}.p = res.p(:);
        out{k,i,2}.mtc = res.mtc;
        out{k,i,2}.lambda = res.lambda;
        
    end
end

