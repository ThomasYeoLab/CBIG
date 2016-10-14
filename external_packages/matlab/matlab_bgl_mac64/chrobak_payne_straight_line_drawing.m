function [X p ei ej] = chrobak_payne_straight_line_drawing(A,varargin)
% CHROBAK_PAYNE_STRAIGHT_LINE_DRAWING Draw planar graphs with straight lines
%
% X = chrobak_payne_straight_line_drawing(A) generates coordinates for each
% vertex such that a planar graph A can be drawn without any edge
% crossings.  This function reports an error if A is not planar.
%
% [X,p,ei,ej] = ... returns additional information.  p is a
% canonical planar ordering of the vertices, and [ei ej] are additional
% edges required to make A a maximal_planar graph.
%
% ... = chrobak_payne_straight_line_drawing(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.is_maximal: is A already a maximal planar graph [{0} | 1]
%
% Note: Be careful with is_maximal=1 and nocheck=1.  If the graph is not
% maximal, then the call will crash Matlab.
%
% Example:
%   [A,xy] = grid_graph(6,5);
%   X = chrobak_payne_straight_line_drawing(A);
%   gplot(A,X,'.-'); hold on; gplot(A,xy*20,'r.-'); hold off
%   % it's still planar, but not obviously a grid!

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2007-10-06: Initial coding
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
            'The graph was not a maximal planar but is_maximal was set.'); end
    end
end

[ei ej p X] = planar_drawing_mex(A,options.is_maximal,0);

