function [d pred f]=astar_search(A,s,h,varargin)
% ASTAR_SEARCH Perform a heuristically guided (A*) search on the graph.
%
% [d pred rank]=astar_search(A,s,h,optionsu) returns the distance map,
% search tree and f-value of each node in an astar_search. 
% The search begins at vertex s.  The heuristic h guides the search, 
% h(v) should be small close to a goal and large far from a goal.  The
% heuristic h can either be a vector with an entry for each vertex in the
% graph or a function which maps vertices to values.
%
% This method works on non-negatively weighted directed graphs.
% The runtime is O((E+V)log(V)).
%
% ... = astar_search(A,u,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.visitor: a visitor to use with the A* search (see Note)
%   options.inf: the value to use for unreachable vertices 
%       [double > 0 | {Inf}]
%   options.target: a special vertex that will stop the search when hit
%       [{'none'} | any vertex number besides the u]; this is ignored if
%       a visitor is set.
%   options.edge_weight: a double array over the edges with an edge
%       weight for each node, see EDGE_INDEX and EXAMPLES/REWEIGHTED_GRAPHS
%       for information on how to use this option correctly
%       [{'matrix'} | length(nnz(A)) double vector]
%   
% Note: You can specify a visitor for this algorithm.  The visitor has the
% following optional functions.
%    vis.initialize_vertex(u)
%    vis.discover_vertex(u)
%    vis.examine_vertex(u)
%    vis.examine_edge(ei,u,v)
%    vis.edge_relaxed(ei,u,v)
%    vis.edge_not_relaxed(ei,u,v)
%    vis.black_target(ei,u,v)
%    vis.finish_vertex(u)
% Each visitor parameter should be a function pointer, which returns 0
% if the search should stop.  (If the function does not return anything, 
% the algorithm continues.)
%
% Example:
%   load graphs/bgl_cities.mat
%   goal = 11; % Binghamton
%   start = 9; % Buffalo
%   % Use the euclidean distance to the goal as the heuristic
%   h = @(u) norm(xy(u,:) - xy(goal,:));
%   % Setup a routine to stop when we find the goal
%   ev = @(u) (u ~= goal);
%   [d pred f] = astar_search(A, start, h, ...
%       struct('visitor', struct('examine_vertex', ev)));

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History 
%  2007-04-20: Added edge weight option
%  2007-07-12: Fixed edge_weight documentation.
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('inf', Inf, 'edge_weight', 'matrix', 'target', 'none');
options = merge_options(options,varargin{:});

edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
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

if check, check_matlab_bgl(A,struct()); end
if trans, A = A'; end

function hi=vec2func(u)
    hi = h(u);
end

if isa(h,'function_handle')
    hfunc = h;
else
    hfunc = @vec2func;
end

if isfield(options,'visitor')
    [d pred f] = astar_search_mex(A,s,target,hfunc,options.inf,edge_weight_opt,options.visitor);
else
    [d pred f] = astar_search_mex(A,s,target,h,options.inf,edge_weight_opt);
end


% end the main function
end
