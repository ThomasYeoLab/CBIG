%% Summary
[G,X] = grid_graph(6,5);
K5 = clique_graph(5);
[is_planar ksubgraph EI]= boyer_myrvold_planarity_test(G);
X = chrobak_payne_straight_line_drawing(G);
test_planar_graph(G);
K = kuratowski_subgraph(K5); % K5 has a Kuratowski subgraph
A = make_maximal_planar(G);
A = make_connected(G);
A = make_biconnected_planar(G);
p = planar_canonical_ordering(G);
[p ei ej] = planar_canonical_ordering(G);
is_kuratowski_graph(K);
is_straight_line_drawing(G,X*20); % round to positions on a grid

%% small empty graphs
for i=0:3
    A = sparse(i,i);
    boyer_myrvold_planarity_test(A);
    X= chrobak_payne_straight_line_drawing(A);
    M = make_maximal_planar(A);
    M = make_connected(A);
    M = make_biconnected_planar(A);
    is_kuratowski_graph(A);
end

%% boyer_myrvold_planarity_test
if boyer_myrvold_planarity_test(clique_graph(5))
    error('matlab_bgl:test_planar','K_5 is not planar');
end
if ~boyer_myrvold_planarity_test(clique_graph(4))
    error('matlab_bgl:test_planar','K_4 is planar');
end
A = clique_graph(5);
A(4,5) = 0;
A(5,4) = 0;
if ~boyer_myrvold_planarity_test(A)
    error('matlab_bgl:test_planar','K_5 - (4,5) is planar');
end

A = clique_graph(5);
[is_planar K] = boyer_myrvold_planarity_test(A);
assert(is_kuratowski_graph(K)==1);
assert(~is_planar,'matlab_bgl:testFailed','K_5 is not planar');
assert(isequal(K,A),'matlab_bgl:testFailed','kuratowski of K_5 = A');

A = clique_graph([3,3]);
[is_planar K] = boyer_myrvold_planarity_test(A);
assert(is_kuratowski_graph(K)==1);
assert(~is_planar,'matlab_bgl:testFailed','K_3,3 is not planar');
assert(isequal(K,A),'matlab_bgl:testFailed','kuratowski of K_3,3 = A');

A = grid_graph(6,5);
[is_planar K] = boyer_myrvold_planarity_test(A);
assert(is_planar==1,'matlab_bgl:testFailed','M(6,5) is planar');
assert(isempty(K),'matlab_bgl:testFailed','no kuratowski for M(6,5)');

% Try some interesting cases
ntrials = 50; nxy=50; 
for i=1:ntrials, 
    x = rand(nxy,1); y=rand(nxy,1);
    tri = delaunay(x,y);
    S = sparse(tri(:,[1 1 2 2 3 3]),tri(:,[2 3 1 3 1 2]),1,nxy, nxy);    
    % S is planar
    [is_planar K EI]= boyer_myrvold_planarity_test(S);
    assert(is_planar==1,'matlab_bgl:testFailed', 'Delaunay graphs are planar');
    n= size(S,1); si= zeros(size(EI.edge_order)); 
    for k=1:n, si(EI.vp(k):EI.vp(k+1)-1)=k; end
    T = sparse(si,EI.edge_order,1);
    assert(isequal(S~=0,T~=0),'matlab_bgl:testFailed', 'Incomplete Embedding Info');
    S(1:5,1:5) = clique_graph(5);
    is_planar= boyer_myrvold_planarity_test(S);    
    assert(is_planar==0,'matlab_bgl:testFailed', 'Delaunay + K5 is not planar');
end

%% is_kuratowski_graph
assert(is_kuratowski_graph(clique_graph(5))==1,'matlab_bgl:testFailed',...
    'K_5 is Kuratowski')
assert(is_kuratowski_graph(clique_graph([3,3]))==1,'matlab_bgl:testFailed',...
    'K_3,3 is Kuratowski')

rand('state',0);
for i=1:10, 
    A=sprand(50,50,5/50); A=A|A'; [isp K]=boyer_myrvold_planarity_test(A);
    if ~is_planar, assert(is_kuratowski_graph(K)==1,'matlab_bgl:testFailed',...
            'invalid Kuratowski graph'); end
end

%% is_straight_line_drawing
X = [0 1; 1 0];
A = clique_graph(2);
assert(is_straight_line_drawing(A,X)==1,'matlab_bgl:testFailed',...
    'the two vertex graph has a straight line drawing');
% is_straight_line_drawing doesn't work right now

A = clique_graph(4);
X = [0 1; 1 0; 0 -1; -1 0];
if is_straight_line_drawing(A,X)==1
    warning('matlab_bgl:bglError',...
    'is_straight_line_drawing failed, old versions of Boost have this problem');
else
    assert(is_straight_line_drawing(A,X)==0,'matlab_bgl:testFailed',...
        'K_4 on a square is not a straight line drawing');
end

%% make_connected
[ei ej] = make_connected(grid_graph(5));
assert(isempty(ei) && isempty(ej), 'matlab_bgl:testFailed',...
    'the grid graph is already connected');

% test 50 random sparse graphs
rand('state',0);
ntrials= 50; n= 50; d= 3;
for i=1:ntrials
    S= sprandsym(n,d/n);
    A= make_connected(S);
    assert(max(components(A))==1, 'matlab_bgl:testFailed',...
        'make_connected did not make the graph connected');
end
    
%% make_biconnected_planar
A = clique_graph(5);
try
    A2 = make_maximal_planar(A);
    assert(1,'matlab_bgl:testFailed','K5 is not planar');
catch end

% test 50 random planar graphs
rand('state',0);
ntrials = 50;
for i=1:ntrials 
    nxy = 50;
    x = rand(nxy,1); y=rand(nxy,1);
    tri = delaunay(x,y);
    S = sparse(tri(:,[1 1 2 2 3 3]),tri(:,[2 3 1 3 1 2]),1,nxy, nxy);
    A = make_biconnected_planar(S);
    assert(nnz(A)>=nnz(S), 'matlab_bgl:testFailed',...
        'make_biconnected_planar removed edges');
    assert(isempty(biconnected_components(A)), 'matlab_bgl:testFailed', ...
        'make_biconnected_planar did not make it biconnected');
end


%% make_maximal_planar
A = clique_graph(5);
try
    A2 = make_maximal_planar(A);
    assert(1,'matlab_bgl:testFailed','K5 is not planar');
catch end

% test 50 random planar graphs
rand('state',0);
ntrials = 50; nxy=50; 
for i=1:ntrials, 
    x = rand(nxy,1); y=rand(nxy,1);
    tri = delaunay(x,y);
    S = sparse(tri(:,[1 1 2 2 3 3]),tri(:,[2 3 1 3 1 2]),1,nxy, nxy);
    A = make_maximal_planar(S);
    assert(nnz(A)>=nnz(S), 'matlab_bgl:testFailed',...
        'make_maximal_planar removed edges');
    [ei ej] = make_maximal_planar(A);
    assert(isempty(ei) && isempty(ej), 'matlab_bgl:testFailed', ...
        'make_maximal_planar added extra edges');
end

%% chrobak_payne_straight_line_drawing
G = grid_graph(10,1);
X = chrobak_payne_straight_line_drawing(G);

rand('state',0);
ntrials = 50; nxy=50; 
for i=1:ntrials, 
    x = rand(nxy,1); y=rand(nxy,1);
    tri = delaunay(x,y);
    S = sparse(tri(:,[1 1 2 2 3 3]),tri(:,[2 3 1 3 1 2]),1,nxy, nxy);
    X= chrobak_payne_straight_line_drawing(S);
    A = make_maximal_planar(S);
    [Y p ei ej]= chrobak_payne_straight_line_drawing(A,'is_maximal',0);    
    assert(isempty(ei) && isempty(ej),'matlab_bgl:testFailed',...
        'is_maximal=0 with maximal graph cannot generate additional edges');
    [Y p ei ej]= chrobak_payne_straight_line_drawing(A,'is_maximal',1);
    assert(isempty(ei) && isempty(ej),'matlab_bgl:testFailed',...
        'is_maximal=1 cannot generate additional edges');
    assert(is_straight_line_drawing(A,X)==1,'matlab_bgl:testFailed',...
        'invalid drawing');
    assert(is_straight_line_drawing(S,X)==1,'matlab_bgl:testFailed',...
        'invalid drawing');
    assert(is_straight_line_drawing(A,Y)==1,'matlab_bgl:testFailed',...
        'invalid drawing');    
    assert(is_straight_line_drawing(S,Y)==1,'matlab_bgl:testFailed',...
        'invalid drawing');
end


%% planar_canonical_ordering(G)
G = grid_graph(6,5);
p = planar_canonical_ordering(G);


%%
% [is_planar ksubgraph EI]= boyer_mvrvold_planarity_test(G);
% X = chrobak_payne_straight_line_drawing(G);
% is_planar_graph(G);
% K = kuratowski_subgraph(G);
% TE = make_maximal_planar(G);
% EE = make_connected(G);
% EE = make_biconnected_planar(G);
% p = planar_canonical_ordering(G);
% [p TE] = planar_canonical_ordering(G);
% is_kurtowski_graph(G);
% is_straight_line_drawing(G,X);