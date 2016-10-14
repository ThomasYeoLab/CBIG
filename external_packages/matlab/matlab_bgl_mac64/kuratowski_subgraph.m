function K = kuratowski_subgraph(A,varargin)
% KURATOWSKI_SUBGRAPH Identify a Kuratowski subgraph 
%
% K = kuratowski_subgraph(A) is just a wrapper around the
% boyer_myrvold_planarity_test to get the subgraph as the first output
% argument.
%
% K is empty if the graph A is planar.  
%
% Example:
%   A = clique_graph([3,3]); % Generate K_3,3
%   K = kuratowski_subgraph(A);
%   isequal(A,K) % K_3,3 is a Kuratowski graph!

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2007-09-29: Initial coding
%%

[is_planar K] = boyer_myrvold_planarity_test(A,varargin{:});
