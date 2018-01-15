function [A xy] = cycle_graph(n,varargin)
% CYCLE_GRAPH Generate the cycle graph of order n
%
% The cycle graph is a simple cycle with n vertices.
%
% [A xy] = cycle_graph(n) generates a cycle graph with n vertices and
% returns the adjacency matrix in A.  The matrix xy stores two-dimensional 
% coordinates for each vertex.
%
% [...] = cycle_graph(n,options) can generate variants on the cycle graph
%   options.directed: return the directed cycle [{0} | 1]
%
% Example:
%   [A xy] = cycle_graph(10);
%   gplot(A,xy);
%
% See also WHEEL_GRAPH, STAR_GRAPH

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-08: Added directed graph option
%  2007-09-29: Changed output to double, fixed for n=0, changed options
%    format
%%

[trans] = get_matlab_bgl_options(varargin{:});
options = struct('directed', '1');
options = merge_options(options,varargin{:});

if n>0
    i = 1:n;
    j = [i(2:end) i(1)];
    A = sparse(i,j,1,n,n);
    if ~options.directed
        A = A|A';
        A = double(A);
    end
    xy = [cos(2*pi*(i./n))' sin(2*pi*(i./n))'];
else
    A = sparse(n,n);
    xy = zeros(0,2);
end
   
if trans && options.directed, A = A'; end
