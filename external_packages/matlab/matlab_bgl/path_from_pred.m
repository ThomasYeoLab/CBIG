function [path P] = path_from_pred(pred,d,varargin)
% PATH_FROM_PRED Convert a predecessor array into a path to a vertex.
%
% path = path_from_pred(pred,d) returns the list of vertices on the path 
% from a source vertex to a vertex d.  The predecessor array is a row 
% vector such that 
%   pred(i) = 0 if vertex i is a source vertex u, 
%   pred(i) = 0 if vertex i has no predecessor associated with it, or
%   pred(i) = j where j preceeds vertex i on the shortest path from u to i.
% The vertex d is a number from 1 to n, where n is the length of the 
% predecessor array.
% 
% The returned path is a 1 by k+1 array where the path from the source to d
% has k edges path lists the order of visting the vertices.  In the case
% that the vertex is the source or unreachable from the source, the path
% just has the single starting vertex.  This reflects one potential problem
% with this function, the source is not uniquely encoded in the predecessor
% array and the source and unreachable vertices are treated equally.
%
% [path P] = path_from_pred(pred,d) also returns the adjacency matrix of the
% directed graph corresponding to the path.  Building this matrix takes
% additional time.
%
% Example:
%    load('graphs/bfs_example.mat');
%    [d dt pred] = bfs(A,1,struct('target', 3));
%    path = path_from_pred(pred,3); % sequence of vertices to vertex 3
%
% See also BFS, DFS, SHORTEST_PATHS

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-04-17: Initial coding
%%

% TODO handle the case of a matrix of predecessors from all_shortest_paths

path = path_from_pred_mex(pred,d);

if nargout == 2
    n = length(pred);
    P = sparse(path(1:end-1), path(2:end), 1, n, n);
end


