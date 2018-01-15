function [p ei ej]=planar_canonical_ordering(A,varargin)
% PLANAR_CANONICAL_ORDERING Compute a planar canonical ordering
%
% See
% http://www.boost.org/doc/libs/1_36_0/libs/graph/doc/planar_canonical_ordering.html
% for a description of a planar canonical ordering
%
% p = planar_canonical_ordering(A) returns the vertex indices in the planar
% canonical ordering.  
%
% [p ei ej] = ... returns the extra edges required to make A a maximal
% planar graph
%
% ... = planar_canonical_ordering(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.is_maximal: is A already a maximal planar graph [{0} | 1]
%
% Note: Be careful with is_maximal=1 and nocheck=1.  If the graph is not
% maximal, then the call will crash Matlab.
%
% Example:
%   p = planar_canonical_ordering(grid_graph(6,5));
%   p(1) % ahh, we draw the first vertex first!

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-10-05: Initial version
%%


[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('is_maximal',0);
options = merge_options(options,varargin{:});
if check
    check_matlab_bgl(A,struct('sym',1)); 
    if options.is_maximal,
        [i j] = make_maximal_planar(A);
        if ~isempty(i), error('matlab_bgl:checkFailed',...
            'The graph was not a maximal planar by is_maximal was set.'); end
    end
end

[ei ej p] = planar_drawing_mex(A,options.is_maximal,1);