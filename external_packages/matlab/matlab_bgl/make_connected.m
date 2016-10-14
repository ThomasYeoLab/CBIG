function varargout= make_connected(A,varargin)
% MAKE_CONNECTED Add edges to construct a connected graph
%
% C = make_connected(A) generates a new graph C with additional edges
% to make C connected.  
%
% [ei,ej] = make_connected(A) returns the additional edges instead.
%
% ... = make_connected(A) takes a set of
% key-value pairs or an options structure.  See set_matlab_bgl_options
% for the standard options. 
%   No additional options for this function
%
% Example:
%   G = sparse(2,2); % empty 2 node graph with 2 components
%   C = make_connected(G)
%   max(components(C))
%   

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

[i j]= planar_edges_mex(A,1,0,0);

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