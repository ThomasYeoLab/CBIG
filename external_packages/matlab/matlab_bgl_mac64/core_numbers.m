function [cn rt] = core_numbers(A,varargin)
% CORE_NUMBERS Compute the core numbers of the vertices in the graph.
%
% [cn rt] = core_numbers(A) returns the core number (cn) of each vertex.
% The core number is the largest integer c such 
% that vertex v exists in a graph where all vertices have degree >= c.
% The vector rt returns the removal time for each vertex.  That
% is, vertex vi was removed at step rt[vi].
%
% This method works on directed graphs but gives the in-degree core number.
% To get the out-degree core numbers, call core_numbers(A').
%
% The runtime is O(E) for unweighted graphs and O((N+M) log N) for weighted
% graphs.  The default is the *unweighted version* which ignores egde
% weights.  For weighted graphs, the definition is the same, but the
% in-degree is the weighted in-degree and is the sum of weight of 
% incoming edges.
%
% To get the out-degree core_numbers, call core_numbers(A') instead;
%
% ... = core_numbers(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.unweighted: an optional switch to perform the weighted 
%       computation [0 | {1}]  
%   options.edge_weight: a double array over the edges with an edge
%       weight for each node, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly
%       [{'matrix'} | length(nnz(A)) double vector]
%
% Note: The default setting for this function is the unweighted computation
% which does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.  This veers from the MatlabBGL
% default of using weighted computations to preserve the standard
% definition of cores for directed and undirected graphs.  
%
% Example: 
%    load graphs/cores_example.mat
%    core_numbers(A)

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-10: Initial version
%  2007-07-11: Updated for weighted cores
%  2007-07-22: Fixed accum array error on Matlab 7.0
%    Added check to not compute size if it isn't used.
%  2007-07-30: Removed size option from output
%    Added removal time output
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('unweighted', 1, 'edge_weight', 'matrix');
options = merge_options(options,varargin{:});

% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix.
edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weights = 1;
    edge_weight_opt = options.edge_weight;
end

if (check)
    % check the values
    if options.unweighted ~= 1 && edge_weights ~= 1
        check_matlab_bgl(A,struct('values',1));
    else
        check_matlab_bgl(A,struct());
    end
end

if trans, A = A'; end

weight_arg = options.unweighted;
if ~weight_arg
    weight_arg = edge_weight_opt;
else
    weight_arg = 0;
end

[cn rt] = core_numbers_mex(A,weight_arg);
% if nargin > 2 && options.unweighted
%     sizes = accumarray([cn+1 ones(length(cn),1)],1);
% end


