function [d pred] = dijkstra_sp(A,u,varargin)
% DIJKSTRA_SP Compute the weighted single source shortest path problem.
%
% Dijkstra's algorithm for the single source shortest path problem only
% works on graphs without negative edge weights.
%
% This method works on weighted directed graphs without negative edge
% weights.
% The runtime is O(V log (V)).
%
% See the shortest_paths function for calling information.  This function 
% just calls shortest_paths(...,struct('algname','dijkstra'));
%
% The options structure can contain a visitor for the Dijkstra algorithm.  
%
% See http://www.boost.org/libs/graph/doc/DijkstraVisitor.html for a 
% description of the events.
% 
% visitor is a struct with the following optional fields
%    vis.initialize_vertex(u)
%    vis.discover_vertex(u)
%    vis.examine_vertex(u)
%    vis.examine_edge(ei,u,v)
%    vis.edge_relaxed(ei,u,v)
%    vis.edge_not_relaxed(ei,u,v)
%    vis.finish_vertex(u)
% Each visitor parameter should be a function pointer, which returns 0
% if the shortest path search should stop.  (If the function does not 
% return anything, the algorithm continues.)
%
% Example:
%    load graphs/clr-25-2.mat
%    dijkstra_sp(A,1)
%
% See also SHORTEST_PATHS, BELLMAN_FORD_SP.

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2008-10-07: Changed options parsing
%%

algname = 'dijkstra';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

[d pred] = shortest_paths(A,u,options);




