load ../graphs/dfs_example.mat
[d dt ft pred] = dfs(A,2);
[ignore order] = sort(dt);
labels(order)


