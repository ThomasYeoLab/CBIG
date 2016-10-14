msgid='matlab_bgl:test_trivial';

%% Test functions on empty or trivial input
% The goal of these tests is to make sure we get error messages or 
% reasonable output (i.e. it doesn't crash) on somewhat malformed or
% trivial input.

% test functions on empty matrices
try
    d = bfs(sparse([]),0);
    error(msgid, 'bfs did not report error');
catch 
end

try
    d = dfs(sparse([]),0);
    error(msgid, 'dfs did not report error');
catch 
end

try
    d = astar_search(sparse([]),0,@(x) x);
    error(msgid, 'astar_search did not report error');
catch 
end

try
    d = shortest_paths(sparse([]), 0);
    error(msgid, 'shortest_paths did not report error');
catch 
end

try
    d = bellman_ford_sp(sparse([]), 0);
    error(msgid, 'bellman_ford_sp did not report error');
catch 
end

try
    d = dag_sp(sparse([]), 0);
    error(msgid, 'dag_sp did not report error');
catch 
end

try
    f = max_flow(sparse([]),0,0);
    error(msgid, 'max_flow did not report error');
catch 
end

try
    p = dominator_tree(sparse([]),0);
    error(msgid, 'lengauer_tarjan_dominator_tree did not report error');
catch 
end
    
D = johnson_all_sp(sparse([]));
D = all_shortest_paths(sparse([]));
D = floyd_warshall_all_sp(sparse([]));
T = mst(sparse([]));
T = kruskal_mst(sparse([]));
T = prim_mst(sparse([]));
cc = components(sparse([]));
bcs = biconnected_components(sparse([]));
c = betweenness_centrality(sparse([]));
c = clustering_coefficients(sparse([]));
ei = edge_weight_index(sparse([]));
m = matching(sparse([]));
m = core_numbers(sparse([]));
v = edge_weight_vector(sparse([]),sparse([]));

X = circle_graph_layout(sparse([]));
X = random_graph_layout(sparse([]));
X = kamada_kawai_spring_layout(sparse([]));
X = fruchterman_reingold_force_directed_layout(sparse([]));
X = gursoy_atun_layout(sparse([]));

is_planar = boyer_myrvold_planarity_test(sparse([]));
[ei ej] = make_connected(sparse([]));
[ei ej] = make_biconnected_planar(sparse([]));
[ei ej] = make_maximal_planar(sparse([]));
p = planar_canonical_ordering(sparse([]));
X = chrobak_payne_straight_line_drawing(sparse([]));
is_kura = is_kuratowski_graph(sparse([]));
try
    is_sldraw = is_straight_line_drawing(sparse([]),zeros(0,2),...
        'nocheck',1,'fix_negative',0);
    error(msgid,'is_straight_line_drawing did not report empty error');
catch 
end