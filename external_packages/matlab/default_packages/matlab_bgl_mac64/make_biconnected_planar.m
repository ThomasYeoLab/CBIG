function varargout= make_biconnected_planar(A,varargin)
% MAKE_BICONNECTED_PLANAR Add edges to construct a biconnected planar graph
%
% M = make_biconnected_planar(A) generates a new graph M with additional edges
% to make M biconnected.  A biconnected planar graph stays connected even
% if any edge is removed.
%
% [ei,ej] = make_biconnected_planar(A) returns the additional edges instead.
%
% ... = make_biconnected_planar(A) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   No additional options for this function
%
% Example:
%   G = grid_graph(6,1); % generate a line graph
%   B = make_biconnected_planar(G);
%   biconnected_components(B)
%   max(components(B))
%   B(1,2) = 0; B(2,1) = 0;
%   max(components(B)) % still one component

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-10-05: Initial version
%%


[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('');
options = merge_options(options,varargin{:});
if check
    check_matlab_bgl(A,struct('sym',1)); 
end

[i j]= planar_edges_mex(A,1,1,0);

if nargout < 2
    if ~isempty(i), 
        [ai aj] = find(A);
        varargout{1} = sparse([ai; i; j], [aj; j; i], 1, size(A,1), size(A,2));
    else
        varargout{1}= A;
    end
else
    varargout{1}= i;
    varargout{2}= j;
end