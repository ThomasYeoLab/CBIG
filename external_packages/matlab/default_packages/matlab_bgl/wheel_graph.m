function [A xy] = wheel_graph(n) 
% WHEEL_GRAPH Construct a wheel graph of order n
%
% The wheel graph is a cycle graph of order n-1 along with an additional
% vertex that connects all the remaining vertices.  (Run the example and it
% will be extremely clear if you are still confused.)
%
% [A xy] = wheel_graph(n) returns the adjacency matrix for the wheel graph
% of order n.  The matrix xy stores two-dimensional coordinates for each 
% vertex.
%
% Example:
%   [A xy] = wheel_graph(10);
%   gplot(A,xy);
%
% See also CYCLE_GRAPH, STAR_GRAPH

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-09-29: Changed output to double, fixed output for i=0,1
%%

if n>1
    i = 1:(n-1);
    j = [i(2:end) i(1)];
    i = [i i];
    j = [j n*ones(1,n-1)];
    A = sparse(i,j,1,n,n);
    A = A|A';
    A = double(A);
else
    i=[];
    A = sparse(n,n);
end

i = i(1:end/2);
if n>0
    xy = [ [cos(2*pi*(i./(n-1))) 0]' [sin(2*pi*(i./(n-1))) 0]' ];
else
    xy = zeros(0,2);
end
