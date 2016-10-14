function [varargout] = prim_mst(A,varargin)
% PRIM_MST Compute a minimum spanning with Prims's algorithm.
%
% Prim's MST algorithm computes a minimum spanning tree for a graph.
%
% This method works on weighted symmetric graphs without negative edge
% weights.
% The runtime is O(E log (V)).
%
% In Boost 1.36.0 and prior versions, Prim's algorithm has a problem with
% incorrect output and diagonal entries.  
%
% See MST for calling information.  This function just calls
% mst(...,struct('algname','prim'));
%   options.root: specify the root vertex for prim's algorithm
%       [{'none'} | any vertex number]
%   options.fix_diag: remove any diagonal entries to get correct output
%       from Prim's algorithm [0 | {1}]; beware this option with the
%       edge_weight option too.
%
% Example:
%    load graphs/clr-24-1.mat
%    prim_mst(A)
%    T = prim_mst(A,struct('root',5));
%
% See also MST, KRUSKAL_MST.

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-23: Initial version
%  2007-04-23: Fixed bug with output without parameters
%  2007-07-08: Fixed documentation typo
%  2008-10-07: Changed options parsing
%%

algname = 'prim';
if ~isempty(varargin), 
    options = merge_options(struct(),varargin{:}); 
    options.algname= algname;
else options = struct('algname',algname); 
end

varargout = cell(1,max(nargout,1));
[varargout{:}] = mst(A,options);

