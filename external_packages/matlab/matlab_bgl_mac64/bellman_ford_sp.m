function [d pred] = bellman_ford_sp(A,u,varargin)
% BELLMAN_FORD_SP Compute the weighted single source shortest path problem.
%
% The Bellman-Ford algorithm for the single source shortest path problem
% works on graphs with negative edge weights.  
%
% See the shortest_paths function for calling information.  This function 
% just calls shortest_paths(...,struct('algname','bellman_ford'));
%
% This method works on weighted directed graphs with negative edge weights.
% The runtime is O(VE).
%
% The options structure can contain a visitor for the Bellman-Ford 
% algorithm.  
%
% See http://www.boost.org/libs/graph/doc/BellmanFordVisitor.html for a 
% description of the events.
% 
% visitor is a struct with the following optional fields
%    vis.initialize_vertex(u)
%    vis.examine_edge(ei,u,v)
%    vis.edge_relaxed(ei,u,v)
%    vis.edge_not_relaxed(ei,u,v)
%    vis.edge_minimized(ei,u,v)
%    vis.edge_not_minimized(ei,u,v)
% Each visitor parameter should be a function pointer, which returns 0
% if the shortest path search should stop.  (If the function does not 
% return anything, the algorithm continues.)
%
% Example:
%    load graphs/kt-6-23.mat
%    d = bellman_ford_sp(A,1);
%
% See also SHORTEST_PATHS, DIJKSTRA_SP.

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2008-10-07: Changed options parsing
%%

algname = 'bellman_ford';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

[d pred] = shortest_paths(A,u,options);


