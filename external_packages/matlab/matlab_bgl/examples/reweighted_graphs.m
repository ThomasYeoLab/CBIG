%% Reweighted graphs in MatlabBGL
% Matlab sparse matrices only support non-zero values. Because the
% structure of the sparse matrix is used to infer the edges of an
% underlying graph this limitation prevents MatlabBGL from trivially 
% addressing graphs with 0-weight edges.
%
% To fix this problem, codes that accept a weighted graph allow the user to
% specify a vector of edge weights for each edge in the graph using the
% optional 'weights' parameter.  Using the 'weights' parameter correctly
% can be difficult due to issues with how edges are indexed in MatlabBGL.

%% Disclaimer
% The details of this section are complicated.  This means their
% implementation is error-prone.  If you get strange behavior, please let
% me know.

%% "I just want the simpliest solution..."
% *new in version 4.0*
% In this section, we'll see the really-easy but somewhat expensive way of
% reweighting a graph.  I'll run through all the cases detailed below with
% the simple code.  If you just need something to work and don't
% necessarily need to know about all the details, this section is for you!

%%
% Let's compute shortest paths in a cycle graph with only one weighted
% edge.  The simple case requires a structural and weight matrix.

%%
% n will be the total size of the graph, and u and v will be the special
% vertices connected with a weight one edge.

n = 8; % it's just an example, so let's make it small.
u = 1;
v = 2;

%%
% These commands create an undirected cycle graph.  The cycle is 
% ... n <-> 1 <-> 2 <-> ... <-> n-1 <-> n <-> 1 ...
% where the weight on every edge is 0 except for the edge between vertex
% u,v.  Notice that the edge list is already symmetric.  
%
% This setup means that while there is a weight 1 edge between u and v, the
% shortest path, or smallest weight path, is actually the path from u,
% circling every vertex except v and so d(v) should be 0.

E = [1:n 2:n 1; 2:n 1 1:n]';
w = [1 zeros(1,n-1) 1 zeros(1,n-1)]';

A = sparse(E(:,1), E(:,2), w, n, n); % create a weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create a structural sparse matrix

%% 
% The relationship between As and A is that As should have a non-zero value
% for every edge, but the _values_ of As will be ignored and the
% computation will proceed with the _values_ in the corresponding spots in
% the matrix A.

%%
% The wrong way to compute shortest paths.

[d pred] = shortest_paths(A,u);
d(v)

%% 
% The right way to compute shortest paths.

[d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));
d(v)

%%
% That's better, d(v) = 0 as expected.


%% A first attempt
% * Correct for version 3.0 *
% A trivial example graph to illustrate the problem that occurs with 0
% weighted graphs occurs even with a simple cycle.  Suppose that the graph
% corresponding to adjacency matrix A is a symmetric cycle where all edges
% have weight 0 except for one edge between vertices (1,2).

%%
% These commands create an undirected cycle graph.  The cycle is 
% ... n <-> 1 <-> 2 <-> ... <-> n-1 <-> n <-> 1 ...
% where the weight on every edge is 0 except for the edge between vertex
% u,v.  Notice that the edge list is already symmetric.  

E = [1:n 2:n 1; 2:n 1 1:n]';
w = [1 zeros(1,n-1) 1 zeros(1,n-1)]';

A = sparse(E(:,1), E(:,2), w, n, n);

%%
% The shortest weighted path between u and v is then through the vertex 
% n because traversing the cycle in the other direction will use the 
% u,v edge of weight 1.  Let's check this with the shortest_paths function.

[d pred] = shortest_paths(A,u);
d(v)

%%
% That is weird, there is a u-v path of length 0 in the graph!  Let's see
% what path the shortest path algorithm chose.

path_from_pred(pred,v)

%%
% The path it chose was from u to v directly, taking the weight 1 edge.
% Let's look at the sparse matrix.

A

%%
% There are only two edges in the matrix corresponding to our symmetric
% weight 1 edge between u and v.  This happens because Matlab removes all 0
% weight edges from the graph.  

%% A first solution
% The solution to the problem is to use the 'edge_weight' optional
% parameter to the shortest_paths function to give it a set of weights 
% to use for each edge.  

help shortest_paths

%% 
% Well, shortest_paths says to read this document, so you are on the right
% track!  It also has a pointer to the function edge_weight_index.  Let's
% look at that function

help edge_weight_index

%%
% This function claims to help us.  It requires building a structural
% matrix which has a non-zero for each edge in the graph.  Let's do that.

As = sparse(E(:,1), E(:,2), 1, n, n)

%%
% Now the matrix has all of the required edges.  According to the
% edge_weight_index function, it returns both a matrix and an index vector.
% The index vector is a way to permute an intelligently ordered set of edge
% weights to the order that MatlabBGL requires the edge weights.  

[ei Ei] = edge_weight_index(As);

full(Ei)
ei

%%
% Now let's create a new edge weight vector for this graph corresponding to
% all the edges we want.  Each non-zero in the matrix should have an
% associated edge weight.  Most the edge weights in this case are 0, so it
% makes it simple.  

ew = zeros(nnz(As),1);
ew(Ei(u,v)) = 1;
ew(Ei(v,u)) = 1;

[d pred] = shortest_paths(As,u,struct('edge_weight', ew(ei)));

path_from_pred(pred,v)

%%
% Excellent, now the shorest path _avoids_ the edge (u,v) as we would
% expect it.  

%% A simplified solution
% The current example is somewhat tedious because we have to create the
% sparse matrix, then create the edge index matrix, and finally create and  
% edit the edge weight array.  
%
% The indexed_sparse function makes the process easier.  
% 
% Recall that using the sparse function directly generated an incorrect
% graph adjacency matrix.

A = sparse(E(:,1), E(:,2), w, n, n)

%%
% The indexed_sparse function is designed as a replacement for sparse where
% the adjacency matrix must be indexed using the edge_weight_index or
% contains 0 weight edges.

help indexed_sparse

%%
% From the documentation of indexed_sparse, the first two return values are
% the structural sparse matrix (As) and the sparse matrix (A) that sparse
% would have returned.  The final two return values are the edge index list
% that edge_weight_index returns as well as the edge index matrix.

% save the old Ei as an example
old_Ei = Ei;

[As A eil Ei] = indexed_sparse(E(:,1), E(:,2), w, n, n);

fprintf('old_Ei = \n\n');
disp(full(old_Ei));
fprintf('Ei = \n\n');
disp(full(Ei))

%% 
% Note that the edge indices changed between the two calls.  The reason for
% this change is that indexed_sparse generates edge indices based on order
% of E(:,1) and E(:,2).  Consequently, this function is much easier to use
% when you already have a set of weighted edges.
%
% In this case, we don't have to create the ew array again!  (Note that the
% call uses eil instead of ei.)

[d pred] = shortest_paths(As,u,struct('edge_weight', w(eil)));

path_from_pred(pred,v)

%% An undirected solution
% The situation for undirected graphs is more complicated.  The trouble
% with the previous solution is that each directed edge had its own weight
% in the vector w.  For an undirected graph, we really want each undirected
% edge to have a single weight, so the natural length of v would be
% nnz(A)/2 instead of nnz(A).
%
% However, MatlabBGL really treats all problems as directed graphs, so it
% will need a vector w of length nnz(A), but that vector should satisfy the
% requirement w(ei1) = w(ei2) if the edges corresponding to ei1 and ei2 are
% (i,j) and (j,i), respectively.  
%
% Again, the edge_weight_index function provides a solution to this
% problem.  We just have to tell edge_weight_index we have an undirected
% graph.

%%
% Let's start with the same sparse matrix

As = sparse(E(:,1), E(:,2), 1, n, n)

%%
% Here we use the edge_weight_index 

[ei Ei] = edge_weight_index(As,struct('undirected',1));
full(Ei) % look at the matrix

%%
% Let's create the edge weight vector
ew = zeros(nnz(As)/2,1); % only half as many zeros here.
ew(Ei(u,v)) = 1; % and we only needed to set one entry to 0

[d pred] = shortest_paths(As,u,struct('edge_weight', ew(ei)));

path_from_pred(pred,v)

%%
% And we get the same output as before!

%% The undirected simplification
% You can probably guess that the simplification for undirected graphs will
% use the indexed_sparse call again too.  

help indexed_sparse

%% 
% From the documentation, we find that indexed_sparse has an option called
% "undirected" which is set to 0 by default.  

[As A eil Ei] = indexed_sparse(E(:,1), E(:,2), w, n, n, struct('undirected',1));

fprintf('Ei = \n\n');
disp(full(Ei))

%%
% In this case, the indexed_sparse routine only issued edge indices that
% were between 1 and 8, rather than 1 and 16 as in the previous case.  

[d pred] = shortest_paths(As,u,struct('edge_weight', w(eil)));

path_from_pred(pred,v)

%% Summary
% The functions that support reweighted edges as of MatlabBGL 3.0 are 
% shortest_paths, all_shortest_paths, dijkstra_sp, bellman_ford_sp, dag_sp,
% betweenness_centrality, astar_search, johnson_all_sp,
% floyd_warshall_all_sp, mst, kruskal_mst, and prim_mst.  Note that
% max_flow does not support these indices.

%%
% The functions that assist working with the edge indices for the
% edge_weight vector are edge_weight_index, indexed_sparse, and
% edge_weight_vector.
