function [A] = clique_graph(n,varargin)
% CLIQUE_GRAPH Generate the clique graph or bipartite clique graph
%
% The clique graph is a graph with as many edges as possible.
%
% A = clique_graph(n) generates a clique with n vertices and returns the
% adjacency matrix A.
%
% A = clique_graph([m n]) generates a bipartite clique with m vertices on
% one side and n vertices on the other side.
%
% A = cycle_graph(...,options) can generate variants of the clique graph
% [...] = cycle_graph(n,options) can generate variants on the cycle graph
%   options.selfloops: add self loops to the graph [{0} | 1]
%
% Example:
%   A = clique_graph(4);
%   test_planar_graph(A);
%   A = clique_graph(5);
%   test_planar_graph(A);

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2007-10-01: Initial coding
%%

options = struct('selfloops', 0);
options = merge_options(options,varargin{:});

if isscalar(n)
    A = ones(n,n);
    if ~options.selfloops,
        A = A - diag(diag(A));
    end
    A = sparse(A);
elseif numel(n) == 2,
    m = n(1);
    n = n(2);
    A = ones(m,n);
    A = spaugment(ones(m,n),0);
    if options.selfloops,
        A = A + speye(m+n);
    end
else 
    error('matlab_bgl:invalidArgument', ...
        'the size option must be a scalar or a pair of numbers');
end
