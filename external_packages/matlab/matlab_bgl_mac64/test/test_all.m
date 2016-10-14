
%
% all_shortest_paths
%
load ../graphs/clr-26-1.mat
floyd_warshall_all_sp(A);

%
% astar_search
% 
load ../graphs/bgl_cities.mat
goal = 11; % Binghamton
start = 9; % Buffalo
% Use the euclidean distance to the goal as the heuristic
h = @(u) norm(xy(u,:) - xy(goal,:));
% Setup a routine to stop when we find the goal
ev = @(u) (u ~= goal);
[d pred f] = astar_search(A, start, h, ...
    struct('visitor', struct('examine_vertex', ev)));

%
% bellman_ford_sp
%
load ../graphs/kt-6-23.mat
d = bellman_ford_sp(A,1);

%
% breadth_first_search
%
load ../graphs/bfs_example.mat
d = test_breadth_first_search(A,1);
d_bfs = bfs(A,1);
if (any(d ~= d_bfs))
   warning('failed breadth_first_search test'); 
end

%
% bfs
%
load ../graphs/bfs_example.mat
d = bfs(A,1);

%
% biconnected components
%
load ../graphs/tarjan-biconn.mat
a = biconnected_components(A);

%
% clustering_coefficients
%
load ../graphs/clique-10.mat
ccfs = clustering_coefficients(A);
if (any(ccfs ~= ones(size(A,1),1)))
    warning('failed clustering_coeffcients test');
end

%
% dag_sp
%
load ../graphs/kt-3-7.mat
d = dag_sp(A,1);
d_bfs = bfs(A,1);
d_bfs(d_bfs < 0) = Inf;
if (any(d ~= d_bfs))
   warning('failed dag_sp test'); 
end



%
% depth_first_search
%
load ../graphs/dfs_example.mat
d = test_depth_first_search(A,1);
d_dfs = dfs(A,1);
if (any(d ~= d_dfs))
   warning('failed depth_first_search test'); 
end


%
% dfs
%
load ../graphs/dfs_example.mat
d = dfs(A,1);
 
%
% dijkstra_sp
%
load ../graphs/clr-25-2.mat
d = dijkstra_sp(A,1);

%
% erdos_reyni
%
A = erdos_reyni(10,0.1);

%
% floyd_warshall_all_sp
%
load ../graphs/clr-26-1.mat
alld = floyd_warshall_all_sp(A);
[alld p] = floyd_warshall_all_sp(A);

%
% johnson_all_sp
%
load ../graphs/clr-26-1.mat
d = johnson_all_sp(A);

%
% kruskal_mst
%
load ../graphs/clr-24-1.mat
T = kruskal_mst(A);

%
% max_flow
%
load ../graphs/max_flow_example.mat
flow = max_flow(A,1,8);
if (flow ~= 4)
    warning('failed max_flow test');
end;

%
% mst
%
load ../graphs/clr-24-1.mat
T = mst(A);


%
% prim_mst
%
load ../graphs/clr-24-1.mat
T = prim_mst(A);

%
% shortest_paths
%
load ../graphs/clr-25-2.mat
[d pred] = shortest_paths(A,1);
[d pred] = shortest_paths(A,1,struct('algname','bellman_ford'));

%
% betweenness_centrality
%

load ../graphs/padgett-florentine.mat
bc = betweenness_centrality(A);
[bc E] = betweenness_centrality(A);
[bc E] = betweenness_centrality(A, struct('unweighted',1));
[bc E] = betweenness_centrality(A, struct('unweighted',1,'ec_list',1));

%
% graph generation
%
A = erdos_reyni(100,.1);
[A xy] = star_graph(10);
[A xy] = cycle_graph(10);
[A xy] = wheel_graph(10);
