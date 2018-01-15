function [A xy] = star_graph(n)
% STAR_GRAPH Generate the star graph of order n
%
% The star graph is a simple star with n vertices.  Vertex n is the center
% vertex.
%
% [A xy] = star_graph(n) generates a star graph with n vertices and
% returns the adjacency matrix in A.  The matrix xy stores two-dimensional 
% coordinates for each vertex.
%
% Example:
%   [A xy] = star_graph(10);
%   gplot(A,xy);
%
% See also WHEEL_GRAPH, CYCLE_GRAPH

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-09-29: Changed output to double
%%


i = 1:n-1;
j = n*ones(1,n-1);
A = sparse(i,j,1,n,n);
A = A|A';
A = double(A);

xy = [[cos(2*pi*(i./(n-1))) 0]' [sin(2*pi*(i./(n-1))) 0]'];
