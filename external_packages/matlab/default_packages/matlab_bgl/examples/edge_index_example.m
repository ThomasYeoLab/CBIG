load ../graphs/bfs_example.mat

[i,j,val] = find(A);

% assign a randon number to each edge in the graph
edge_rand = rand(num_edges(A),1);

Av = sparse(i,j,edge_rand, size(A,1), size(A,2));

ee = @(ei,u,v) fprintf('examine_edge   %2i,   %1i,   %1i,   %4f,   %4f\n', ...
    ei, u, v, edge_rand(ei), Av(u,v));

fprintf('\n-------\n');
fprintf('Example of using the transposed edge index incorrectly\n');
fprintf('-------\n');
fprintf('               ei,   u,   v,     er(ei),     A(u,v)\n');
breadth_first_search(A,1,struct('examine_edge',ee));

% build a tranposed edge index to edge index map
Aind = sparse(i,j,1:num_edges(A),size(A,1), size(A,2));
[i,j,trans_ei_to_ei] = find(Aind');


ee = @(ei,u,v) fprintf('examine_edge   %2i,   %1i,   %1i,   %4f,   %4f\n', ...
    ei, u, v, edge_rand(trans_ei_to_ei(ei)), Av(u,v));

fprintf('\n-------\n');
fprintf('Example of using the transposed edge index correctly,\n');
fprintf('but with a different value specified for (u,v) and (v,u)\n')
fprintf('-------\n');
fprintf('               ei,   u,   v,     er(ei),true er(ei)\n');
breadth_first_search(A,1,struct('examine_edge',ee));



% build a transposed edge index to edge index map respecting how each edge
% is undirected


[i,j,val] = find(triu(A,1));

edge_rand = rand(num_edges(A)/2,1);

Av = sparse([i; j], [j; i], [edge_rand; edge_rand], size(A,1), size(A,2));

% build a tranposed edge index to edge index map
Aind = sparse([i; j],[j; i],[1:num_edges(A)/2 1:num_edges(A)/2], size(A,1), size(A,2));
[i,j,trans_ei_to_ei] = find(Aind');





ee = @(ei,u,v) fprintf('examine_edge   %2i,   %1i,   %1i,   %4f,   %4f\n', ...
    ei, u, v, edge_rand(trans_ei_to_ei(ei)), edge_rand(Aind(u,v)));

fprintf('\n-------\n');
fprintf('Example of using the transposed edge index correctly,\n');
fprintf('but a single value for both (u,v) and (v,u) edges\n')
fprintf('-------\n');
fprintf('               ei,   u,   v,     er(ei),true er(ei)\n');
breadth_first_search(A,1,struct('examine_edge',ee));


