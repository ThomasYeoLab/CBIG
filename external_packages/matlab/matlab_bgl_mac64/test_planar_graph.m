function is_planar = test_planar_graph(G)
% TEST_PLANAR_GRAPH Test if a graph is planar
%
% This function is just a convinence wrapper around
% boyer_myrvold_planarity_test.
%
% Example:
%   A = clique_graph(5);
%   test_planar_graph(A)
%   A(4,5)=0; A(5,4)=0; % make K_5 into a planar graph
%   test_planar_graph(A)

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2007-09-29: Initial coding
%%

is_planar = boyer_myrvold_planarity_test(G);
