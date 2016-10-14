function v=edge_weight_vector(As,A,varargin)
% EDGE_WEIGHT_VECTOR Returns input for the edge_weight option for a weighted
% matrix
%
% Given a structural and weighted matrix pair (As,A), this function returns
% input for the edge_weight option.  The structural matrix As can have
% arbitrary non-zero values, but the non-zero structure of A must be a
% subset of the non-zero structure of As.  In terms of graphs, think about
% it as: the weights come from A, but the edges come from As.
%
% Example:
%   n = 8; u = 1; v = 2;
%   E = [1:n 2:n 1; 2:n 1 1:n]';
%   w = [1 zeros(1,n-1) 1 zeros(1,n-1)]';
%   A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
%   As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
%   [d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));
%   d(v)

% David Gleich
% Copyright, Stanford University, 2008

%% History
%  2008-09-24: Initial coding
%%
[trans check] = get_matlab_bgl_options(varargin{:});
if check
    check_matlab_bgl(A,struct()); 
    if any(size(As)~=size(A)), error('matlab_bgl:edge_weight_vector', ...
            'The structral and weight matrix must be the same size'); end
end
n = size(A,1);
if trans, [j i] = find(As');
else [i j] = find(As);
end
Ai = sparse(i,j,1:nnz(As),n,n);
inds = nonzeros(Ai.*(A~=0));
if length(inds) ~= nnz(A), error('matlab_bgl:edge_weight_vector', ...
    'The matrix A cannot have additional non-zeros not in As.'); end
v = zeros(nnz(As),1);
v(inds) = nonzeros(A);
if check && ~isequal(A,sparse(j,i,v,n,n)')
    error('matlab_bgl:edge_weight_vector','error in weight function');
end

%% Todo
% Check that is actually works with the is trans option, I'm not convinced
% the test coverage works.