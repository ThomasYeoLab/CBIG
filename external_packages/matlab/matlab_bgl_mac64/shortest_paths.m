function [d pred] = shortest_paths(A,u,varargin)
% SHORTEST_PATHS Compute the weighted single source shortest path problem.
%
% [d pred] = shortest_paths(A,u) returns the distance (d) and the predecessor
% (pred) for each of the vertices along the shortest path from u to every
% other vertex in the graph.  
% 
% ... = shortest_paths(A,u,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.algname: the algorithm to use 
%       [{'auto'} | 'dijkstra' | 'bellman_ford' | 'dag']
%   options.inf: the value to use for unreachable vertices 
%       [double > 0 | {Inf}]
%   options.target: a special vertex that will stop the search when hit
%       [{'none'} | any vertex number besides the u]; target is ignored if
%       visitor is set.
%   options.visitor: a structure with visitor callbacks.  This option only
%       applies to dijkstra or bellman_ford algorithms.  See dijkstra_sp or
%       bellman_ford_sp for details on the visitors.
%   options.edge_weight: a double array over the edges with an edge
%       weight for each edge, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly
%       [{'matrix'} | length(nnz(A)) double vector]
%
% Note: if you need to compute shortest paths with 0 weight edges, you must
% use an edge_weight vector, see the examples for details.
%
% Note: 'auto' cannot be used with 'nocheck' = 1.  The 'auto' algorithm
% checks if the graph has negative edges and uses bellman_ford in that
% case, otherwise, it uses 'dijkstra'.  In the future, it may check if the
% graph is a dag and use 'dag'.  
%
% Example:
%    load graphs/clr-25-2.mat
%    shortest_paths(A,1)
%    shortest_paths(A,1,struct('algname','bellman_ford'))
%
% See also DIJKSTRA_SP, BELLMAN_FORD_SP, DAG_SP

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-19: Initial coding
%  2007-04-18: Added edge_weight option.
%  2007-04-19: Added target option.
%    Added additional error checks.
%  2007-07-12: Fixed edge_weight documentation
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('algname', 'auto', 'inf', Inf, 'edge_weight', 'matrix', ...
    'target', 'none');
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

if strcmp(options.target,'none')
    target = 0; % a flag used to denote "no target" to the mex
elseif isa(options.target, 'double')
    target = options.target;
else
    error('matlab_bgl:invalidParameter', ...
        'options.target is not ''none'' or a vertex number.');
end

if check
    % check the values of the matrix
    check_matlab_bgl(A,struct('values',edge_weights ~= 1));
    
    if edge_weights && nnz(A) ~= length(edge_weight_opt)
        error('matlab_bgl:invalidParameter', 'the vector of edge weights must have length nnz(A)');
    end
    
    % set the algname
    if (strcmpi(options.algname, 'auto'))
        if edge_weights
            mv = min(edge_weights);
        else
            mv = min(min(A));
        end
        
        if (mv < 0)
            options.algname = 'bellman_ford';
        else
            options.algname = 'dijkstra';
        end
    else
        % check the data provided to match the algorithm
        if strcmpi(options.algname, 'dijkstra')
            if edge_weights
                mv = min(edge_weight_opt);
            else
                mv = min(min(A));
            end
            if mv < 0
                error('matlab_bgl:invalidParameter', ...
                    'dijkstra''s algorithm cannot be used with negative edge weights.');
            end
        end
    end
    
else
    if (strcmpi(options.algname, 'auto'))
        error('shortest_paths:invalidParameter', ...
            'algname auto is not compatible with no check');       
    end
end

if options.inf < 0, error('options.inf must be larger than 0'); end

if trans, A = A'; end

if isfield(options,'visitor')
    [d pred] = matlab_bgl_sp_mex(A,u,target,lower(options.algname),options.inf,...
        edge_weight_opt, options.visitor);
else
    [d pred] = matlab_bgl_sp_mex(A,u,target,lower(options.algname),options.inf,...
        edge_weight_opt);
end

