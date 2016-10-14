function [l c] = path_histogram(G,varargin)
% PATH_HISTOGRAM Compute a histogram of all shortest paths in graph G
%
% [l c] = path_histogram(G) computes all shortest paths in A one at a time
% and forms the histogram of the shortest distances between all vertices.  
%
% [l c] = path_histogram(G,struct('sample',N)) uses N random samples instead of an
% exact solution.  (N=0 uses exact solution)
%
% Options:
%   options.weighted: use Dijkstra algorithm for weighted paths [{0} | 1]
%   options.intcount: use integer histogram bins [0 | {1}]
%   options.nbins: number of bins when intcount=0 [{50} | positive integer]
%   options.sample: use N iid samples instead of all sources
%     [{0} | positive integer]
%   options.verbose: output status for long runs [{0} | 1]
%
% Example:
%  load('graphs/cs-stanford.mat');
%  [l c] = path_histogram(A);
%  bar(l,c,[min(l),max(l)],'hist');

%
% David Gleich
% Copyright, Stanford University, 2008
%

% History
% 2008-03-14: Initial coding by David

options=struct('weighted',0,'intcount',1,'nbins',50','sample',0,'istrans',0,...
    'verbose',0);
if ~isempty(varargin), options = merge_structs(varargin{1}, options); end

if ~options.istrans, G = G'; end
bfsoptions=struct('istrans',1);

n=size(G,1); verb=options.verbose; 

N = options.sample;
if N==0, vs=1:n; N = n;
else vs=ceil(n*rand(options.sample,1)); 
end

% sample from all vertices
H=sparse(n,1); t0=clock; p=1; Np=min(N,100);
for vi=1:N
    if verb && vi==(p*floor(N/Np) + max(mod(N,Np)-(Np-p),0))
        dt=etime(clock,t0); p=p+1;
        fprintf('  %5.1f : %9i of %9i; etr=%7f sec\n', ...
            100*(p-1)/Np, vi, N, N*dt/vi-dt); 
    end
    di=bfs(G,vi,bfsoptions);
    di=di(di>0);             % only use reachable points
    H=H+accumarray(di,1,[n 1],@sum,0,true);
    bfsoptions.nocheck=1;    % don't repeat options check
end
l = find(H); c = nonzeros(H);

end

% internal copy of merge_structs
function S=merge_structs(A,B)
S = A; fn = fieldnames(B);
for ii = 1:length(fn), if (~isfield(A, fn{ii})), S.(fn{ii}) = B.(fn{ii}); end, end
end
