function [D,P] = floyd_warshall_all_sp(A,varargin)
% FLOYD_WARSHALL_ALL_SP Compute the weighted all-pairs shortest path problem.
%
% The Floyd-Warshall algorithm for the all-pairs shortest path problem 
% works only on graphs without negative edge weights.  This method should
% be used over the Johnson algorithm for dense graphs.  
%
% This algorithm can return the predecessor matrix.
%
% This method works on weighted directed graphs.
% The runtime is O(V^3).
%
% See the shortest_paths function for calling information.  This function 
% just calls all_shortest_paths(...,struct('algname','floyd_warshall'));
%
% Example:
%    load graphs/clr-26-1.mat
%    floyd_warshall_all_sp(A)
%
% See also ALL_SHORTEST_PATHS, JOHNSON_ALL_SP.

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2008-04-02: Added documenation for predecessor matrix
%  2008-10-07: Changed options parsing
%%

algname = 'floyd_warshall';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

if nargout > 1, [D,P] = all_shortest_paths(A,options);
else D = all_shortest_paths(A,options);
end
