function C = approx_multiway_cut(A,vs)
% APPROX_MULTIWAY_CUT Solve a 2-approximation to the multi-way cut problem
% 
% C = approx_multiway_cut(A,vs)
%
% Outputs C, the set of edges cut in a 2-approximation to the multiway cut
% problem.  The multiway-cut problem is to find a minimum cost set of edges
% to disconnect all the vertices in vs from each other.
%
% The non-zero values contain the weight of each edge.
%
% The input A must be a symmetric graph.

if (~isequal(A,A'))
    error('approx_multiway_cut:invalidParameter',...
        'the matrix must be symmetric.');
end;

if (min(min(A)) < 0)
    error('approx_multiway_cut:invalidParameter',...
        'the matrix cannot contain negative weights.');
end;

n = size(A,1);

% this should be larger than any conceivable flow...
int_infinity = sum(sum(A))+2*sum(sum(A(vs,:)))+1;

% initial the cut to nothing.
C = sparse(n,n);

% Get A as an edge list...
[i j] = find(A);

for kk=1:length(vs)
    v = vs(kk);
    others = setdiff(vs,v);
    
    % Each flow problem add a fake sink as the n+1 vertex
    Aflow = A;
    Aflow(others,n+1) = int_infinity*ones(length(others),1);
    Aflow(n+1,:) = sparse(n+1,1);
    
    % solve the max-flow problem
    [flow ci] = max_flow(Aflow,v,n+1);
    
    % remove the last (fake) entry from the cut.
    ci = ci(1:end-1);
    
    % construct a value over the edges that is 0 except on the cut, we know
    % all values are positive, so just take the absolute value
    vc = abs(v.*(ci(i)-ci(j)))./2;
    
    % add the set of edges to the cut by constructing a sparse matrix with
    % only the cut edges.
    C = C+ sparse(i, j, vc, n,n);
end;