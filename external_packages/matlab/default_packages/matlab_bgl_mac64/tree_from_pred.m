function T = tree_from_pred(pred)
% TREE_FROM_PRED Convert a predecessor array into a tree 
%
% T = tree_from_pred(pred) returns the predecessor array as an
% unweighted adjacency matrix
%   pred(i) = 0 if vertex i is a source vertex u, 
%   pred(i) = 0 if vertex i has no predecessor associated with it, or
%   pred(i) = j where j is the parent of vertex i in the tree.
% 
% Example:
%    load('graphs/dominator_tree_example.mat');
%    p = lengauer_tarjan_dominator_tree(A,1);
%    T = tree_from_pred(p);
%
% See also BFS, DFS, LENGAUER_TARJAN_DOMINATOR_TREE

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-12: Initial version
%  2008-09-23: Reformatted common section. 
%%

n = length(pred);
mask = pred~=0;
I = 1:n; 
T = sparse(I(mask),pred(mask),1,n,n);
