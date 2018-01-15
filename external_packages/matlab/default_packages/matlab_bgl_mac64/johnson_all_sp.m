function D = johnson_all_sp(A,varargin)
% JOHNSON_ALL_SP Compute the weighted all-pairs shortest path problem.
%
% Johnson's algorithm for the all-pairs shortest path problem 
% works only on graphs without negative edge weights.  This method should
% be used over the Floyd-Warshall algorithm for sparse graphs.  
%
% This method works on weighted directed graphs.
% The runtime is O(VE log(V)).
%
% See the shortest_paths function for calling information.  This function 
% just calls all_shortest_paths(...,struct('algname','johnson'));
%
% Example:
%    load graphs/clr-26-1.mat
%    johnson_all_sp(A)
%
% See also ALL_SHORTEST_PATHS, FLOYD_WARSHALL_ALL_SP.


% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2008-09-24: Code cleanup
%  2008-10-07: Changed options parsing
%%

algname = 'johnson';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

D = all_shortest_paths(A,options);




