function n=num_edges(A)
% NUM_EDGES The number of edges in a graph.
%
% n = num_edges(A) returns the number of edges in graph A.  
%
% For symmetric/undirected graphs, the number of edges returned is twice 
% the number of undirected edges.
%
% Example:
%    load graphs/dfs_example.mat
%    n = num_edges(A)
%
% See also NUM_VERTICES

% David Gleich
% Copyright, Stanford University, 2006-2008

n = nnz(A);
