function D = dijkstra_all_sp(G,optionsu)
% DIJKSTRA_ALL_SP Compute all pairs shortest path using repeated Dijkstras
%
% This algorithm runs a series of Dijkstra's shortest path problems to
% compute the all pairs shortest path matrix. 
%
% Using one of the other algorithms is preferable, however, this algorithm
% is useful as a base when there are memory concerns with the other
% algorithms.  For example, you could easily parallelize this algorithm to
% compute a subset of shortests paths on each processor and then aggregate
% the results.
%
% D = dijkstra_all_sp(G) produces identical output to
% floyd_warshall_all_sp(G).  
%
% Example:
%    load graphs/clr-26-1.mat
%    D1 = floyd_warshall_all_sp(A)
%    D2 = dijkstra_all_sp(A)
%
% See also ALL_SHORTEST_PATHS, JOHNSON_ALL_SP, FLOYD_WARSHALL_ALL_SP.

%
% TODO: Make the code work for dijkstra or bellman_ford
%

%
% 13 July 2007
% Implement the algorithm correct for transposed input
%

% First, transpose the graph so that Dijkstra's won't have to at each
% iteration.
if exist('optionsu','var')
    options = optionsu;
else
    options = struct('istrans',0);
end

if ~options.istrans
    G = G';
end

% Next, check to make sure a Dijkstra's call works so we don't have to
% check on the inner loop.
options.istrans=1;
d_i = dijkstra_sp(G,1,options);

options.nocheck=1;
D = zeros(size(G));
for vi=2:size(G,1)
    D(:,vi-1) = d_i;
    d_i = dijkstra_sp(G,vi,options);
end
% save the last one
D(:,end) = d_i;

if ~options.istrans
    D = D';
end


