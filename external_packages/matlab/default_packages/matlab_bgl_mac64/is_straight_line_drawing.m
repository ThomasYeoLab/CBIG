function is_sldrawing = is_straight_line_drawing(A,X,varargin)
% IS_STRAIGHT_LINE_DRAWING Test if coordinates are a straight line drawing
%
% is_sldrawing = is_straight_line_drawing(A,X) determines if the
% coordinates of each vertex in X allow the graph A to be drawn with
% straight lines and no edges crossings.  is_sldrawing=1 if this is true
% and =0 otherwise.
%
% Internally, the coordinates X cannot be negative and must be integers,
% the default options automatically adjust the coordinates.
%
% This function uses a bucket-sort and may take large amounts of memory if
% the coordinates in X are not well-placed.
%
% ... = is_straight_line_drawing(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   options.fix_negative: automatically adjust negative layouts [0 | {1}]
%
% Example:
%   X = [0 1; 1 0];
%   is_straight_line_drawing(clique_graph(2),X)
%   X = [0 1; 1 0; 0 -1; -1 0];
%   is_straight_line_drawing(clique_graph(4),X)
%   % Oops, some version of the BGL have an error than doesn't correctly
%   % detect the second case.



% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-10-05: Initial version
%%


[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('fix_negative',1);
options = merge_options(options,varargin{:});
if check
    check_matlab_bgl(A,struct('sym',1)); 
    floorerr = max(max(abs(floor(X)-X)));
    if floorerr>eps(1), warning('matlab_bgl:roundWarning',...
      'after floor, the values in X changed by %g, possible incorrect output',...
      floorerr);
    end
end

if check || options.fix_negative
    minX = min(min(X));
    maxX = max(max(X));
    if check && maxX-minX>10*size(A,1) && maxX-minX>1e7,
        warning('matlab_bgl:memoryWarning',...
            'is_straight_line_drawing uses a bucket sort that may run out of memory');
    end
    if options.fix_negative && minX < 0
        X = X-minX;
    end
end

is_sldrawing = planar_test_mex(A,2,X);