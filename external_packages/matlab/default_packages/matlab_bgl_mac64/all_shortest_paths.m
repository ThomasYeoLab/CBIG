function [D,P] = all_shortest_paths(A,varargin)
% all_shortest_paths Compute the weighted all pairs shortest path problem.
%
% D = all_shortest_paths(A) returns the distance matrix D for all vertices
% where D(i,j) indicates the shortest path distance between vertex i and
% vertex j.  
%
% For the Floyd-Warshall algorithm, this function can return the 
% predecessor matrix as well:
%   [D P]=all_shortest_paths(A,struct('algname','floyd_warshall'));
% returns P so that the P(i,j) is the node preceeding j on the path from
% i to j.  To build the path between (i,j), use the commands
%   p=[]; while j~=0, p(end+1)=j; j=P(i,j); end; p=fliplr(p);
% 
% ... = all_shortest_paths(A,u,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.algname: the algorithm to use 
%       [{'auto'} | 'johnson' | 'floyd_warshall']
%   options.inf: the value to use for unreachable vertices 
%       [double > 0 | {Inf}]
%   options.edge_weight: a double array over the edges with an edge
%       weight for each node, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly
%       [{'matrix'} | length(nnz(A)) double vector]
%
% Note: 'auto' cannot be used with 'nocheck' = 1.  The 'auto' algorithms
% checks the number of edges in A and if the graph is more than 10% dense,
% it uses the Floyd-Warshall algorithm instead of Johnson's algorithm.
%
% Example:
%    load graphs/clr-26-1.mat
%    all_shortest_paths(A)
%    all_shortest_paths(A,struct('algname','johnson'))
%
% See also JOHNSON_ALL_SP, FLOYD_WARSHALL_ALL_SP.

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-19: Initial version
%  2006-05-31: Added full2sparse check
%  2007-03-01: Added option for predecessor matrix from floyd_warshall
%  2007-04-20: Added edge weight option
%  2007-07-08: Fixed typos in strings and documentation
%    Removed fixes for the Johnson algorithm
%  2007-07-12: Fixed edge_weight documentation.
%  2007-07-21: Fixed divide by 0 error in check for algorithm type
%  2008-04-02: Added documenation for predecessor matrix
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('algname', 'auto', 'inf', Inf, 'edge_weight', 'matrix');
options = merge_options(options,varargin{:});

% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix.
%edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weight_opt = options.edge_weight;
end

if check
    % check the values of the matrix
    check_matlab_bgl(A,struct('values',1));
    
    % set the algname
    if strcmpi(options.algname, 'auto')
        nz = nnz(A);
        if (nz/(numel(A)+1) > .1)
            options.algname = 'floyd_warshall';
        else
            options.algname = 'johnson';
        end
    end
else
    if strcmpi(options.algname, 'auto')
        error('all_shortest_paths:invalidParameter', ...
            'algname auto is not compatible with no check');       
    end
end

if trans, A = A'; end

if nargout > 1
    [D,P] = matlab_bgl_all_sp_mex(A,lower(options.algname),options.inf,edge_weight_opt);
    P = P';
else
    D = matlab_bgl_all_sp_mex(A,lower(options.algname),options.inf,edge_weight_opt);
end

if trans, D = D'; end

