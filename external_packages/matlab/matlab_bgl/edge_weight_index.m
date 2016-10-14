function [ei Ei] = edge_weight_index(A,varargin)
% EDGE_WEIGHT_INDEX Build a conformal matrix of edge index values for a graph.
%
% [eil Ei] = edge_weight_index(As) returns a vector where 
%   As(i,j) not= 0 implies Ei(i,j) not= 0 and Ei(i,j) = eil(i)
% for an integer value of eil(i) that corresponds to the edge index value
% passed in the visitors.  
%
% The input matrix A should be a structural matrix with a non-zero value
% for each edge.  The matrix Ei gives an index for each edge in the graph,
% and the vector eil will reorder a vector of edge weights to an appropriate
% input for 'edge_weight' parameter of a function call.
%
% The edge_weight_index function assists writing codes that use the
% edge_weight parameter to reweight a graph based on a vector of weights
% for the graph or using the ei parameter from an edge visitor.  It is 
% critical to obtain high performance when
%
% i) constructing algorithms that use 0 weighted edges
% ii) constructing algorithms that change edge weights often.
%
% See the examples reweighted_edges and edge_index_visitor for more
% information.
%
% ... = edge_weight_index(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%    options.undirected: output edge indices for an undirected graph [{0} | 1]
%        see Note 1.
%
% Note 1: For an undirected graph, the edge indices of the edge corresponding
% to (u,v) and (v,u) are the same.  Consequently, Ei is a symmetric matrix,
% using this option allows only one value for an undirected edge.
%
% Example:
%   load('graphs/bfs_example.mat');
%   [eil Ei] = edge_weight_index(A,struct('undirected',1));
%   edge_rand = rand(num_edges(A)/2,1);
%   [iu ju] = find(triu(A,0));
%   Av = sparse(iu,ju,edge_rand,size(A,1),size(A,1)); Av = Av + Av';
%   ee = @(ei,u,v) fprintf('examine_edge %2i, %1i, %1i, %4f, %4f, %4f\n', ...
%               ei, u, v, edge_rand(eil(ei)), Av(u,v), edge_rand(Ei(u,v)));
%   breadth_first_search(A,1,struct('examine_edge', ee));
%
% See also INDEXED_SPARSE

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-13: Changed input options to use undirected as the option name.
%  2007-07-24: Fixed example
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('undirected', 0);
options = merge_options(options, varargin{:});
if check, check_matlab_bgl(A,struct('sym',options.undirected == 0)); end

if options.undirected
    % compute the numer of edges
    [i,j] = find(triu(A,0));
    ne = length(i);
    
    diag_mask = i~=j;
    eis = 1:ne;
    
    Ei = sparse([i;j(diag_mask)],[j;i(diag_mask)],[eis eis(diag_mask)], size(A,1), size(A,2));
else
    [i,j] = find(A);
    ne = length(i);
    
    Ei = sparse(i,j,1:ne,size(A,1), size(A,2));
end

if trans
    ei = nonzeros(Ei');
else
    ei = nonzeros(Ei);
end

