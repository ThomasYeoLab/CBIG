function varargout= make_maximal_planar(A,varargin)
% MAKE_MAXIMAL_PLANAR Add edges to construct a maximal planar graph
%
% M = make_maximal_planar(A) generates a new graph M with additional edges
% so that adding any other edge will make a non-planar graph.  
%
% [ei,ej] = make_maximal_planar(A) returns the additional edges instead.
%
% ... = make_maximal_planar(A) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   No additional options for this function
%
% Example:
%   G = grid_graph(6,5);
%   M = make_maximal_planar(G);
%   test_planar_graph(M) % it's planar!
%   M(9,20) = 1; M(20,9) = 1; 
%   test_planar_graph(M) % it's not planar


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

[i j]= planar_edges_mex(A,1,1,1);

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