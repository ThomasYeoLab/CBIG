load ../graphs/max_flow_example.mat
max_flow(A,1,8)
[flow cut R F] = max_flow(A,1,8);
full(R)     


