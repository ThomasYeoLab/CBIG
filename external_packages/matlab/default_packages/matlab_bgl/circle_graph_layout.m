function X = circle_graph_layout(G,radius)
% CIRCLE_LAYOUT Layout the vertices of a graph on a circle
% 
% X = circle_layout(G) generates a layout of graph with vertices uniformly
% placed on the circle.  This function does no interesting layout and just
% places the vertices around the circle in order.
%
% X = circle_layout(G,radius) places the vertices with a radius other than
% 1.0.
%
% Example:
%   G = cycle_graph(6);
%   X = circle_graph_layout(G);
%   gplot(G,X);

% David F. Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-09-25: Initial coding
%%

if ~exist('radius','var') || isempty(radius), radius = 1.0; end;
pi = 3.14159;
n = num_vertices(G);
X = zeros(n,2);
X(:,1) = radius*cos( (0:n-1)'*2*pi/n );
X(:,2) = radius*sin( (0:n-1)'*2*pi/n );
