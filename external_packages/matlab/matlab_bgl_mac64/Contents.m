% Graph Algorithms
%
% Searching through a graph
% bfs                       - Breadth first search
% dfs                       - Depth first search
% astar_search              - Heuristic astar graph search
% breadth_first_search      - Breadth first search with visitors
% depth_first_search        - Depth first search with visitors
%
% Shortest Path Algorithms
% shortest_paths            - Single source shortest path wrapper
% all_shortest_paths        - All pairs shortest path wrapper
% dijkstra_sp               - Dijkstra's shorest path algorithm
% bellman_ford_sp           - Bellman-Ford shortest path algorithm
% dag_sp                    - Shortest path on directed acyclic graph
% johnson_all_sp            - Johnson all pairs shortest path algorithm
% floyd_warshall_all_sp     - Floyd-Warshall all pairs shortest path alg
%
% Minimum Spanning Tree
% mst                       - Minimum spanning tree wrapper
% kruskal_mst               - Kruskal's minimum spanning tree algorithm
% prim_mst                  - Prim's minimum spanning tree algorithm
%
% Connected Components
% components                - Connected components of a graph
% biconnected_components    - Biconnected connected components of a graph
%
% Flow Algorithms
% max_flow                  - Solve a maximum flow problem
% edmunds_karp_max_flow     - Edmunds-Karp max flow algorithm
% kolmogorov_max_flow       - Kolmogorov's max flow algorithm
% push_relabel_max_flow     - Goldberg's push-relabel max flow algorithm
%
% Layouts
% circle_graph_layout       - Simple layout of vertices on a circle
% random_graph_layout       - Random layout of vertices in plane or lattice
% kamada_kawai_spring_layout- Spring based graph layout
% gursoy_atun_layout        - Topology based graph layout
% fruchterman_reingold_force_directed_layout - Force directed graph layout
%
% Matchings
% matching                  - Compute a maximum cardinality matching
% edmonds_maximum_cardinality_matching - Edmonds' algorithm for matching 
% maximal_matching          - Compute maximal (not maximum) matchings
% test_matching             - Test if a matching is maximum cardinality
%
% Statistics
% betweenness_centrality    - Betweeness centrality scores for all nodes
% clustering_coefficients   - Clustering coefficients for all nodes
% core_numbers              - Compute in-degree core numbers for all nodes
% lengauer_tarjan_dominator_tree - Compute a dominator tree for a graph
% num_edges                 - The number of edges in a graph
% num_vertices              - The number of vertices in a graph
% topological_order         - Compute a topological order for a dag
% test_dag                  - Test if a graph is directed and acyclic
%
% Graphs 
% clique_graph              - Generates a clique or bipartite clique
% cycle_graph               - Generates a cycle graph
% erdos_reyni               - Generates an erdos_reyni, or Gnp, graph
% grid_graph                - Generate a grid or hypergrid graph
% star_graph                - Generates a star graph
% wheel_graph               - Generates a wheel graph 
%
% Visitors
% combine_visitors          - Produce a new combination visitor
%
% Utilities
% edge_weight_index         - Convert between graphs and edge indices
% edge_weight_vector        - Generate edge_weight vectors from matrices
% indexed_sparse            - Generate a sparse matrix with edge indices
% path_from_pred            - Convert a predecessor array to a path
% tree_from_pred            - Convert predecessor array to a tree
% 
% Examples and Demos
% EXAMPLES/RED_BLACK        - Compute a red-black ordering of a matrix
% examples/record_alg       - Use visitors to show how an algorithm works
% examples/reweighted_graphs    - Show how reweighted graphs work
% examples/core_numbers_example - Demonstrate core numbers
% examples/planar_graph     - A few planar graph examples
% examples/new_in_3_0       - New features in version 3.0 
% examples/new_in_4_0       - New features in version 4.0
% examples/layouts          - Simple demonstrations of layout algorithms
%
% Options
% set_matlab_bgl_default    - Set default options


% future functions...
% sample_paths              - Computes path statistics

% Matrix Orderings
% graph_perm                 - Graph permutation wrapper
% reverse_cuthill_mckee_perm - Reverse Cuthill-McKee ordering
% minimum_degree_perm        - Minimum degree ordering
% king_perm                  - King ordering
% sloan_perm                 - Sloan ordering

% David Gleich
% Copyright, Stanford University, 2006-2008
