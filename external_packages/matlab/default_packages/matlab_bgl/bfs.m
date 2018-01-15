function [d dt pred] = bfs(A,u,varargin)
% BFS Compute the breadth first search order.
%
% [d dt pred] = bfs(A,u) returns the distance to each vertex (d) and the  
% discover time (dt) in a breadth first search starting from vertex u.
%    d(i) = dt(i) = -1 if vertex i is not reachable from vertex u.
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
%
% This method works on directed graphs.
% The runtime is O(V+E).
%
% ... = bfs(A,u,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.target: a special vertex that will stop the search when hit
%       [{'none'} | any vertex number besides the u]
%
% Note: this function does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.
%
% Example:
%    load graphs/bfs_example.mat
%    d = bfs(A,1)
%
% See also DFS

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

options = struct('target', 'none');
options = merge_options(options,varargin{:});

if check, check_matlab_bgl(A,struct()); end

if strcmp(options.target,'none')
    target = 0; % a flag used to denote "no target" to the mex
elseif isa(options.target, 'double')
    target = options.target;
else
    error('matlab_bgl:invalidParameter', ...
        'options.target is not ''none'' or a vertex number.');
end

if (trans) 
    A = A'; 
end

[d dt pred] = bfs_mex(A,u,target);


