function [As,A,eil,Ei] = indexed_sparse(i,j,v,m,n,varargin)
% INDEXED_SPARSE Create a sparse matrix with indexed edges.
%
% [As,A,eil,Ei] = indexed_sparse(i,j,v,m,n) creates a sparse matrix A just
% like A = sparse(i,j,v,m,n).  However, indexed_sparse returns additional
% information.  The matrix As is a structural matrix for A which
% corresponds to As = sparse(i,j,1,m,n).  Thus, As(i,j) != 0 for all edges.
% The vector eil is a permutation for the vector v, such that v(eil) is the
% correct input for the edge_weight parameter.  The matrix Ei lists the
% index of each edge in the vector v, so that 
% A = sparse(j,i,v(nonzeros(Ei)),m,n)' unless options.istrans = 0.  
%
% This function handles the case when v(k) == 0.  For v(k) = 0,
% A(i(k),j(k)) = 0, but As(i(k),j(k)) = 1, and the vector v(eil) provides
% an appropriate input to the edge_weight parameter for all the algorithms.
%  
% See the examples reweighted_edges for more information.
%
% ... = indexed_sparse(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%    options.undirected: output edge indices for an undirected graph [{0} | 1]
%      See the note about undirected inputs.
%
% Note (Undirected inputs): If options.undirected = 1, the input to the
% graph still must contain both undirected edges and the corresponding
% weight.  
%
% Example:
%   % see example/reweighted_edges 

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-13: Changed input options to use undirected as the option name.
%  2008-10-07: Changed options parsing
%%

[trans check]  = get_matlab_bgl_options(varargin{:});

options = struct('undirected', 0);
options = merge_options(options, varargin{:});

symmetric = options.undirected;

Ei = accumarray([i j], 1:length(v), [m n], @min, 0, true);
if symmetric, Ei = min(Ei,Ei'); end

A = sparse(i, j, v, m, n);
As = spones(Ei);

if trans, eil = nonzeros(Ei'); 
else eil = nonzeros(Ei);
end

if check
    % make sure there were no elements in i,j repeated 
    % that would cause problems
    if length(v) ~= nnz(Ei)
        warning('matlab_bgl:indexed_sparse', ...
            'duplicated elements detected, using the first occurrence as the index');
    end
    
    if symmetric
        [ei ej] = find(As');
        if ~isequal(A,sparse(ei,ej,v(eil),m,n))
            warning('matlab_bgl:indexed_sparse', ...
                'the output matrix failed a symmetry test.  This indicates a non-symmetric value array v.');
        end
    end
end


