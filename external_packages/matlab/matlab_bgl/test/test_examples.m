function test_examples

%% Code examples
msgid='matlab_bgl:test_examples';

% all_shortest_paths
load('../graphs/clr-26-1.mat');
all_shortest_paths(A);
all_shortest_paths(A,struct('algname','johnson'));

[D P]=all_shortest_paths(A,struct('algname','floyd_warshall'));
i=4; j=1; p=[]; while j~=0, p(end+1)=j; j=P(i,j); end; p=fliplr(p);

% astar_search
load('../graphs/bgl_cities.mat');
goal = 11; % Binghamton
start = 9; % Buffalo
% Use the euclidean distance to the goal as the heuristic
h = @(u) norm(xy(u,:) - xy(goal,:));
% Setup a routine to stop when we find the goal
ev = @(u) (u ~= goal);
[d pred f] = astar_search(A, start, h, ...
    struct('visitor', struct('examine_vertex', ev)));

% bellman_ford_sp
load('../graphs/kt-6-23.mat');
d = bellman_ford_sp(A,1);

% betweenness_centrality
load('../graphs/padgett-florentine.mat');
betweenness_centrality(A);

% bfs
load('../graphs/bfs_example.mat');
d = bfs(A,1);

% biconnected_components
load('../graphs/tarjan-biconn.mat');
biconnected_components(A);

% breadth_first_search
% see (dist_uv_bfs below)    
load('../graphs/bfs_example.mat');
d2 = dist_uv_bfs(A,1,3);

% clustering_coefficients
load('../graphs/clique-10.mat');
clustering_coefficients(A);

% combine_visitors
vis1 = struct();
vis1.examine_vertex = @(u) fprintf('vis1: examine_vertex(%i)\n', u);
vis2 = struct();
vis2.examine_vertex = @(u) fprintf('vis2: examine_vertex(%i)\n', u);
combined_vis = combine_visitors(vis1, vis2);
load('../graphs/bfs_example.mat');
breadth_first_search(A,1,combined_vis);

% components
load('../graphs/dfs_example.mat');
components(A);
     
% core_numbers
load('../graphs/cores_example.mat');
cn = core_numbers(A);

% cycle_graph
[A xy] = cycle_graph(10);
gplot(A,xy);

% dag_sp
load('../graphs/kt-3-7.mat');
dag_sp(A,1);

% depth_first_search
load('../graphs/dfs_example.mat');
dist_uv_dfs(A,1,3);

% dfs
load('../graphs/dfs_example.mat');
d = dfs(A,1);

% dijkstra_sp
load('../graphs/clr-25-2.mat')
[d pred] = dijkstra_sp(A,1);

% edge_weight_index
load('../graphs/bfs_example.mat');
[eil Ei] = edge_weight_index(A,struct('undirected',1));
edge_rand = rand(num_edges(A)/2,1);
[iu ju] = find(triu(A,0));
Av = sparse(iu,ju,edge_rand,size(A,1),size(A,1)); Av = Av + Av';
ee = @(ei,u,v) fprintf('examine_edge %2i, %1i, %1i, %4f, %4f, %4f\n', ...
            ei, u, v, edge_rand(eil(ei)), Av(u,v), edge_rand(Ei(u,v)));
breadth_first_search(A,1,struct('examine_edge', ee));

% edge weight vector
n = 8; u = 1; v = 2;
E = [1:n 2:n 1; 2:n 1 1:n]';
w = [1 zeros(1,n-1) 1 zeros(1,n-1)]';
A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
[d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));

% edmonds_maximum_cardinality_matching
load('../graphs/matching_example.mat');
m=edmonds_maximum_cardinality_matching(A) ;
sum(m>0)/2;                % matching cardinality, should be 8

% edmunds_karp_max_flow
load('../graphs/max_flow_example.mat');
f=edmunds_karp_max_flow(A,1,8);
    
% erdos_reyni
A = erdos_reyni(100,0.05);

% floyd_warshall_all_sp
load('../graphs/clr-26-1.mat');
D=floyd_warshall_all_sp(A);

% grid_graph
[A xy] = grid_graph(5,10);
gplot(A,xy);
A = grid_graph(2*ones(1,10)); % compute 10d hypercube

    
% indexed_sparse

% johnson_all_sp
load('../graphs/clr-26-1.mat');
D=johnson_all_sp(A);

% kolmogorov_max_flow
load('../graphs/max_flow_example.mat');
kolmogorov_max_flow(A,1,8);

% kruskal_mst
load('../graphs/clr-24-1.mat');
T=kruskal_mst(A);
     
% lengauer_tarjan_dominator_tree
load('../graphs/dominator_tree_example.mat');
p = lengauer_tarjan_dominator_tree(A,1);

% matching
load('../graphs/matching_example.mat');
[m,v] = matching(A);
[m,v] = matching(A,struct('augmenting_path','none'));

% matching
load('../graphs/matching_example.mat');
m = matching(A); 
sum(m>0)/2;                % matching cardinality, should be 8
[m,v] = matching(A,struct('augmenting_path','none')); % not maximum matching

% max_flow
load('../graphs/max_flow_example.mat');
f=max_flow(A,1,8);

% maximal_matching
load('../graphs/matching_example.mat');
m=maximal_matching(A);
sum(m>0)/2;                % maximal matching cardinality, should be < 8
mmax=matching(A);
sum(mmax>0)/2;             % maximum matching cardinality, should be 8 

% mst
load('../graphs/clr-24-1.mat');
T=mst(A);

% num_edges
load('../graphs/dfs_example.mat');
n = num_edges(A);

% num_vertices
A = sparse(ones(5));
n = num_vertices(A);

% path_from_pred
load('../graphs/bfs_example.mat');
[d dt pred] = bfs(A,1,struct('target', 3));
path = path_from_pred(pred,3); % sequence of vertices to vertex 3

% prim_mst
load('../graphs/clr-24-1.mat');
prim_mst(A);
T = prim_mst(A,struct('root',5)); % root the tree at vertex e

% push_relabel_max_flow
load('../graphs/max_flow_example.mat');
f=push_relabel_max_flow(A,1,8);


% shoretst_paths
load('../graphs/clr-25-2.mat');
shortest_paths(A,1);
shortest_paths(A,1,struct('algname','bellman_ford'));

% star_graph
[A xy] = star_graph(10);
gplot(A,xy);

% test_dag
n = 10; A = sparse(1:n-1, 2:n, 1, n, n); % construct a simple dag
test_dag(A);
A(10,1) = 1; % complete the cycle
test_dag(A);

% test_matching
load('../graphs/matching_example.mat');
[m_not_max,v] = matching(A,struct('augmenting_path','none'));
test_matching(A,m_not_max);

% toplogical_order
load('../graphs/bfs_example.mat');
d = bfs(A,1);

% tree_from_pred
load('../graphs/dominator_tree_example.mat');
p = lengauer_tarjan_dominator_tree(A,1);
T = tree_from_pred(p);

% wheel graph
[A xy] = wheel_graph(10);
gplot(A,xy);
n = 10;
A = cycle_graph(n);
[d dt ft pred] = dfs(A,1,struct('target',3));

%% Graph layout algorithms
msgid = 'test_examples:layout';

% circle_graph_layout
G = cycle_graph(6);
X = circle_graph_layout(G);
gplot(G,X);

% fruchterman_reingold_force_directed_layout
G = grid_graph(6,5);
X = fruchterman_reingold_force_directed_layout(G);
gplot(G,X);

% gursoy_atun_layout
G1 = cycle_graph(5000,struct('directed',0));
X1 = gursoy_atun_layout(G1,'topology','heart');
G2 = grid_graph(50,50);
X2 = gursoy_atun_layout(G2,'topology','square');
G3 = grid_graph(50,50);
X3 = gursoy_atun_layout(G3,'topology','circle');
subplot(1,3,1); gplot(G1,X1,'k'); subplot(1,3,2); gplot(G2,X2,'k');
subplot(1,3,3); gplot(G3,X3,'k');

% kamada_kawai_spring_layout
G = grid_graph(6,5);
X = kamada_kawai_spring_layout(G);
gplot(G,X);

% random_graph_layout
G = cycle_graph(1500);
X = random_graph_layout(G);
gplot(G,X); hold on; plot(X(:,1),X(:,2),'r.'); hold off;
% Layout on the grid
X = random_graph_layout(G,int32([0 0 5 5])); % random grid layout
gplot(G,X); grid on; hold on; plot(X(:,1),X(:,2),'r.'); hold off;

%% Planar graph algorithms
msgid = 'test_examples:planar';

% boyer_myrvold_planarity_test
G = grid_graph(6,5);
assert(boyer_myrvold_planarity_test(G)==1,'msgid','Grid(6,5) is planar');
K5 = clique_graph(5);
assert(boyer_myrvold_planarity_test(K5)==0,'msgid','K5 is not planar');

% chrobak_payne_straight_line_drawing
[A,xy] = grid_graph(6,5);
X = chrobak_payne_straight_line_drawing(A);
gplot(A,X,'.-'); hold on; gplot(A,xy*20,'r.-'); hold off
% it's still planar, but not obviously a grid!

% test_planar_graph
A = clique_graph(5);
assert(test_planar_graph(A)==0,msgid,'K5 is not planar');
A(4,5)=0; A(5,4)=0; % make K_5 into a planar graph
assert(test_planar_graph(A)==1,msgid,'K5-edge is planar');

% kuratowski_subgraph
A = clique_graph([3,3]); % Generate K_3,3
K = kuratowski_subgraph(A);
assert(isequal(A,K)==1,msgid,'K_3,3 is a Kuratowski graph!');

% make_maximal_planar
G = grid_graph(6,5);
M = make_maximal_planar(G);
assert(test_planar_graph(M)==1,msgid,'Maximal(Grid(6,5)) is planar');
M(9,20) = 1; M(20,9) = 1; 
assert(test_planar_graph(M)==0,msgid,'Maximal(Grid(6,5)) + edge is not planar');

% make_biconnected_planar
G = grid_graph(6,1); % generate a line graph
B = make_biconnected_planar(G);
assert(isempty(biconnected_components(B)),msgid,'make_biconnected failed');
assert(max(components(B))==1,msgid,'biconnected graphs have one component');
B(1,2) = 0; B(2,1) = 0;
assert(max(components(B))==1,msgid,'biconnected graphs - edge have one component');

% make_connected
G = sparse(2,2); % empty 2 node graph with 2 components
C = make_connected(G);
assert(max(components(C))==1,msgid,'make_connected failed');

% is_straight_line_drawing
X = [0 1; 1 0];
assert(is_straight_line_drawing(clique_graph(2),X)==1,msgid,...
    'straight line test failed');
X = [0 1; 1 0; 0 -1; -1 0];
is_straight_line_drawing(clique_graph(4),X);

% is_kuratowski_graph
assert(is_kuratowski_graph(clique_graph(4))==0,msgid,'K4 is not kuratowski');
assert(is_kuratowski_graph(clique_graph(5))==1,msgid,'K5 is kuratowski');
assert(is_kuratowski_graph(clique_graph([3,3]))==1,msgid,'K_3,3 is kuratowski');

% planar_canonical_ordering
p = planar_canonical_ordering(grid_graph(6,5));
assert(p(1)==1,msgid,'Grid(6,5) can draw vertex 1 first');


end

%% Accessory functions

    function dv=dist_uv_bfs(A,u,v)
      vstar = v;
      dmap = ipdouble(zeros(size(A,1),1));
      function stop=on_tree_edge(ei,u,v)
        dmap(v) = dmap(u)+1;
        stop = (v ~= vstar);
      end
      breadth_first_search(A,u,struct('tree_edge',@on_tree_edge));
      dv = dmap(v);
    end
    
    function dist_uv_dfs(A,u,v)
      vstar = v;
      dmap = ipdouble(zeros(size(A,1),1));
      function stop=on_tree_edge(ei,u,v)
        dmap(v) = dmap(u)+1;
        stop = (v ~= vstar);
      end
      breadth_first_search(A,u,struct('tree_edge',@on_tree_edge));
    end

