function is_kgraph = is_kuratowski_graph(A,varargin)
% IS_KURATOWSKI_GRAPH Test if a graph can be collapsed to K_3,3 or K_5
%
% is_k = is_kuratowski_graph(A) checks if A can be reduced to K_3,3 or K_5
% by repeated edge contractions.  If so, then is_k=1, otherwise is_k=0.
%
% ... = is_straight_line_drawing(A,...) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   No additional options for this function
%
% Example:
%   is_kuratowski_graph(clique_graph(4)) % false, K_4 is not kuratowski
%   is_kuratowski_graph(clique_graph(5)) % true, K_5 is kuratowski
%   is_kuratowski_graph(clique_graph([3,3])) % true, K_3,3 is kuratowski

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-10-05: Initial version
%%


[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct();
options = merge_options(options,varargin{:});

if check, check_matlab_bgl(A,struct('sym',1)); end

is_kgraph = planar_test_mex(A,1);