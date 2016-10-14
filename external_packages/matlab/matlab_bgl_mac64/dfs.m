function [d dt ft pred] = dfs(A,u,varargin)
% DFS Compute the depth first search times.
%
% [d dt ft pred] = dfs(A,u) returns the distance (d), the discover (dt) and
% finish time (ft) for each vertex in the graph in a depth first search 
% starting from vertex u.
%   d = dt(i) = ft(i) = -1 if vertex i is not reachable from u
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
% 
% ... = dfs(A,u,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.full: compute the full dfs instead of the dfs of
%      the current component (see Note 1) [{0} | 1]
%   options.target: a special vertex that will stop the search when hit
%       (see Note 2) [{'none'} | any vertex number besides the u]
%
% Note 1: When computing the full dfs, the vertex u is ignored, vertex 1 is
% always used as the starting vertex.  
%
% Note 2: When target is specified, the finish time array only records the
% finish time for all vertices that actually finished.  The array will then
% be undefined on a significant portion of the graph, but that doesn't
% indicate the vertices are unreachable; they just haven't been reached
% yet.
%
% Note: this function does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.
%
% Example:
%    load graphs/dfs_example.mat
%    d = dfs(A,1)
%
% See also BFS

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-19: Initial version
%  2006-05-31: Added full2sparse check
%  2007-04-19: Added target option
%  2008-10-07: Changed options parsing
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('target', 'none', 'full', 0);
options = merge_options(options, varargin{:});

if check, check_matlab_bgl(A,struct()); end

if strcmp(options.target,'none')
    target = 0; % a flag used to denote "no target" to the mex
elseif isa(options.target, 'double')
    target = options.target;
else
    error('matlab_bgl:invalidParameter', ...
        'options.target is not ''none'' or a vertex number.');
end

if trans, A = A'; end

[d dt ft pred] = dfs_mex(A,u,target,options.full);

