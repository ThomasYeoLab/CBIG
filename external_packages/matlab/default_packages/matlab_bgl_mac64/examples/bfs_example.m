load ../graphs/bfs_example.mat
[d dt pred] = bfs(A,2);
[ignore order] = sort(dt);
labels(order)
treeplot(pred);

