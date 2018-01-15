function [ci sizes] = components(A,varargin)
% COMPONENTS Compute the connected components of a graph.
%
% [ci sizes] = components(A) returns the component index vector (ci) and
% the size of each of the connected components (sizes).  The number of
% connected components is max(components(A)).  The algorithm used computes
% the strongly connected components of A, which are the connected
% components of A if A is undirected (i.e. symmetric).  
%
% This method works on directed graphs.
% The runtime is O(V+E), the algorithm is just depth first search.
%
% ... = components(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   There are no additional options for this function.
%
% Note: this function does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.
%
% Example: 
%    load('graphs/dfs_example.mat');
%    components(A)
%
% See also DMPERM, BICONNECTED_COMPONENTS

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-21: Initial version
%  2006-05-31: Added full2sparse check
%  2006-11-09: Fixed documentation typo.
%  2007-07-08: Code cleanup
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end
if check, check_matlab_bgl(A,struct()); end
if trans, A = A'; end

[ci sizes] = components_mex(A);


