function [is_planar ksubgraph EI]=boyer_myrvold_planarity_test(A,varargin)
% BOYER_MYRVOLD_PLANARITY_TEST Test a graph for planarity
%
% is_planar = boyer_myrvold_planaity_test(A) yields 1 if A is a planar
% graph or 0 otherwise.  A planar graph can be drawn on the xy plane with
% no edge crossings.
%
% [is_planar K] = ... identifies a Kuratowski subgraph of A when A is not
% planar.
% [is_planar K EI] = ... computes the planar embedding edge order in
% EI.edge_order.  EI.vp gives the starting and ending entries in EI.edge_order
% for the order of edges from each vertex.  (EI is the compressed
% sparse row representation of the matrix A with the order of the column
% indices permuted to the planar graph embedding order.)
%
% ... = boyer_myrvold_planaity_test(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   No additional options for this function
%
% Example:
%   G = grid_graph(6,5);
%   boyer_myrvold_planarity_test(G) % G is planar
%   K5 = clique_graph(5);
%   boyer_myrvold_planarity_test(K5) % K5 is not planar

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2007-10-06: Initial coding
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct();
options = merge_options(options,varargin{:});

if nargout <= 1
    is_planar = planar_test_mex(A,0);
else 
    if nargout <= 2
        [is_planar ki kj] = planar_test_mex(A,0);
    else
        [is_planar ki kj eip eie] = planar_test_mex(A,0);
        EI = struct('vp',eip,'edge_order',eie);
    end
    if ~isempty(ki)
        ksubgraph = sparse(ki,kj,1,size(A,1),size(A,2));
    else
        ksubgraph = sparse([]);
    end
    ksubgraph = ksubgraph|ksubgraph'; % placed here to get the type right    
end