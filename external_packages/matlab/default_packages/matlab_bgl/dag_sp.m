function [d pred] = dag_sp(A,u,varargin)
% DAG_SP Compute the weighted single source shortest path problem.
%
% The DAG shortest path algorithm for the single source shortest path
% problem only works on directed acyclic-graphs (DAGs).  
%
% If the graph is not a DAG, the results are undefined.  In the future, the
% function may throw an error if the graph is not a DAG.
%
% See the shortest_paths function for calling information.  This function 
% just calls shortest_paths(...,struct('algname','dag'));
%
% This algorithm works on weighted directed acyclic graphs.
% The runtime is O(V+E)
%
% ... = clustering_coefficients(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%    There are no additional options for this function.
%
% Example:
%    load graphs/kt-3-7.mat
%    dag_sp(A,1)
%
% See also SHORTEST_PATHS

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2008-10-07: Changed options parsing
%%

algname = 'dag';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

[d pred] = shortest_paths(A,u,options);




