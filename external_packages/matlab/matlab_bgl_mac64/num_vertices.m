function n=num_vertices(A)
% NUM_VERTICES The number of vertices in a graph.
%
% n = num_vertices(A) returns the number of vertices in graph A.  
%
% Example:
%    A = sparse(ones(5));
%    n = num_vertices(A);
%
% See also NUM_EDGES

% David Gleich
% Copyright, Stanford University, 2006-2008

n = size(A,1);
