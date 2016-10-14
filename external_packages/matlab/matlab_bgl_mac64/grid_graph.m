function [A coords] = grid_graph(varargin)
% GRID_GRAPH Generate a grid graph or hypergrid graph 
% 
% [A xy] = grid_graph(m,n) generates a grid graph with m vertices along the
% x axis and n vertices along the y axis.  The xy output gives the 2d
% coordinates of each vertex. 
% [A xyz] = grid_graph(m,n,k) generates a grid graph with m vertices along
% the x axis, n vertices along the y axis, and k vertices along the z axis.
% The xyz output gives the 3d coordinates of each vertex.
% [A X] = grid_graph(d1, d2, ..., dk) generates a hypergrid graph in k
% dimensions with d1 vertices along dimension 1, d2 vertices along
% dimension 2, ..., and dk vertices along the final dimension.  The X
% output gives the k dimensional coordinates of each vertex.
% [A X] = grid_graph([d1, d2, ..., dk]) is an alternate input scheme.
%
% Example:
%   [A xy] = grid_graph(5,10);
%   gplot(A,xy);
%   A = grid_graph(2*ones(1,10)); % compute 10d hypercube

% David Gleich
% Copyright, Stanford University, 2007-2008

%% History
%  2007-07-13: Initial version
%%

k = length(varargin);

if k==1 && numel(varargin{1}) > 1
    varargin = num2cell(varargin{1});
    k = length(varargin);
else
    if any(cellfun('prodofsize',varargin)~=1)
        error('matlab_bgl:invalidArgument',...
            'please specific the size of dimension in the arguments');
    end
end

dims=cell(k,1);
dim_size = cell2mat(varargin(end:-1:1));

A = sparse(1);
for ii=1:length(dim_size)
    n2 = dim_size(ii);
    dims{ii} = linspace(0,1,n2);
    
    %A = kron(A,line_graph(n2)) + kron(line_graph(n2),speye(size(A,1)));
    A = kron(speye(size(A,1)), line_graph(n2)) + kron(A,speye(n2));
end

% retarded, but it'll have to do...
if (nargout > 1)
    if (k == 1)
        % this case is easy enough
        n = dim_size(1);
        coords = (0:n-1)'./(n-1);
    else
        % Matlab's ndgrid does something different for one dimension
        % we also have to reverse the output and get each component of the
        % output.  
        coords_cell = cell(k,1);
        cmdstr = '[';
        for ii=1:k
            cmdstr = [cmdstr sprintf('coords_cell{%i} ',k-ii+1)];      %#ok
        end
        cmdstr = [cmdstr '] = ndgrid(dims{:});'];
        eval(cmdstr);

        coords = zeros(size(A,1),k);
        for ii=1:k
            coords(:,ii) = reshape(permute(coords_cell{ii},k:-1:1),size(A,1),1);
        end
    end
end

function [Al] = line_graph(n)
e = ones(n,1);
Al = spdiags([e e], [-1 1], n, n);
