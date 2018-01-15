function [varargout] = kruskal_mst(A,varargin)
% KRUSKAL_MST Compute a minimum spanning with Kruskal's algorithm.
%
% The Kruskal MST algorithm computes a minimum spanning tree for a graph.
%
% This method works on weighted symmetric graphs.
% The runtime is O(E log (E)).
%
% See the mst function for calling information.  This function just calls
% mst(...,struct('algname','kruskal'));
%
% Example:
%    load graphs/clr-24-1.mat
%    kruskal_mst(A)
%
% See also MST, PRIM_MST.

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2008-09-24: Code cleanup
%%

algname = 'kruskal';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

varargout = cell(1,max(nargout,1));

[varargout{:}] = mst(A,options);

