%% Recording algorithm behavior with MatlabBGL
% In this example, we will write a simple visitor that outputs an
% algorithm's behavior.  The algorithm we will examine is dijkstra_sp.  To
% examine the runtime behavior we will use a visitor which outputs a string
% every time a function is called.


%% Setup
% To begin, we load a graph.

load ../graphs/clr-25-2.mat

%%
% Next, let's check the documentation to see which functions to implement 
% for the visitor

help dijkstra_sp

%%
% The help states that dijkstra_sp allows visitors functions for
% initialize_vertex, discover_vertex, examine_vertex, examine_edge,
% edge_relaxed, edge_not_relaxed, and finish_vertex.
%
% Rather than implementing 7 functions ourselves, we define two helper
% functions.  These helper functions return functions themselves.  There is
% one helper that returns a vertex visitor function and one helper than
% returns an edge visitor function.

vertex_vis_print_func = @(str) @(u) ...
    fprintf('%s called on %s\n', str, char(labels{u}));
edge_vis_print_func = @(str) @(ei,u,v) ...
    fprintf('%s called on (%s,%s)\n', str, char(labels{u}), char(labels{v}));

%% 
% These anonymous functions return functions themselves.

ev_func = vertex_vis_print_func('examine_vertex');
ev_func(1)

%% 
% I hope you see how these functions are useful in saving quite a bit of
% typing.

%% Calling dijkstra_sp
% We are almost done.  Now, we just have to setup the visitor structure to
% pass to the dijkstra_sp call.

vis = struct();
vis.initialize_vertex = vertex_vis_print_func('initialize_vertex');
vis.discover_vertex = vertex_vis_print_func('discover_vertex');
vis.examine_vertex = vertex_vis_print_func('examine_vertex');
vis.finish_vertex = vertex_vis_print_func('finish_vertex');
vis.examine_edge = edge_vis_print_func('examine_edge');
vis.edge_relaxed = edge_vis_print_func('edge_relaxed');
vis.edge_not_relaxed = edge_vis_print_func('edge_not_relaxed');

%%
% With the visitor setup, there is hardly any work left.  

dijkstra_sp(A,1,struct('visitor', vis));

%% Understanding the output
% To understand the output, we find it helpful to have a copy of
% Introduction to Algorithms by Cormen, Leiserson, and Rivest.  The source
% for the graph is Figure 25-2 in that book and the authors use the graph
% to illustrate how Dijkstra's algorithm runs.  In particular, Figure 25-5
% shows a sample run of Dijkstra's algorithm.
%
% Perhaps the first thing to notice is that the initialize vertex visitor
% is never called.  This results from an error in the MatlabBGL and Boost
% documentation.  Once it is resolved, we will update the MatlabBGL
% documentation to match the Boost graph library.
%
% The results: discover_vertex is called before examine_vertex.  For the 
% edges, examine_edge is always called before either edge_relaxed
% or edge_not_relaxed.  The edges that are relaxed are the shaded edges in
% Figure 25-5.
%
% Finally, finish vertex is called on a vertex after all of its edges have
% been examined and possibly relaxed.  