function [a,C] = biconnected_components(A,varargin)
% BICONNECTED_COMPONENTS Compute the biconnected components and
% articulation points for a symmetric graph A.
%
% [a C] = biconnected_components(A) returns a list of articulation points 
% a and the component graph C where each non-zero indicates the connected
% component of the edge.  That is, C is a matrix with the same non-zero
% structure as A, but with the values replaced with the index of the
% biconnected component of that edge.  The vector a is a list of
% articulation points in the graph.  Articulation points are vertices that
% belong to more than one biconnected component.  Removing an articulation
% point disconnects the graph.
%
% If C is not requested, it is not built.
%
% This method works on undirected graphs.
% The runtime is O(V+E), the algorithm is just depth first search.
%
% ... = biconnected_components(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   There are no additional options for this function.
%
% Note: the input to this function must be symmetric, so this function
% ignores the 'notrans' default option and never transposes the input.
%
% Note: this function does not depend upon the non-zero values of A, but
% only uses the non-zero structure of A.
%
% Example:
%    load graphs/tarjan-biconn.mat
%    biconnected_components(A)
%
% See also COMPONENTS

% David Gleich
% Copyright, Stanford University, 2006-2008

%% History
%  2006-04-19: Initial version
%  2006-05-31: Added full2sparse check
%%

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

if trans, end

if check
    % make sure the matrix is symmetric
    check_matlab_bgl(A,struct('sym',1));
end;

% the graph has to be symmetric, so trans doesn't matter.

if (nargout > 1)
    [a ci] = biconnected_components_mex(A);
    
    % convert the indices from the graph back into a new matrix.
    [i j] = find(A);
    C = sparse(i,j,ci,size(A,1),size(A,1));

    C = max(C,C');    
else
    a = biconnected_components_mex(A);
end;

% 0 a indicates it isn't an articulation point.
a = a(a > 0);




