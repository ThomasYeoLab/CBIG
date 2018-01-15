function test_statistics

%%
msgid = 'matlab_bgl:test_statistics';

%% betweenness_centrality

n = 10;
A = cycle_graph(n);

% the centrality index of all vertices in a cycle is the same
bc = betweenness_centrality(A,struct('unweighted',1));
if any(bc-bc(1))
    error(msgid, 'betweenness_centrality returned incorrect values for a cycle');
end

% make sure we toss an error when the graph has logical weights
try 
    bc = betweenness_centrality(A);
    error(msgid, 'betweenness_centrality did not report an error');
catch
end

% make sure the edge centrality graphs are the same in a few cases
A = sparse(1:n-1, 2:n, 1, n, n);
[bc,Ec] = betweenness_centrality(A);
if any(any(spones(A) - spones(Ec)))
    error(msgid, 'different non-zero structure in edge centrality matrix');
end

[bc,Ec] = betweenness_centrality(A,struct('istrans',1));
if any(any(spones(A) - spones(Ec)))
    error(msgid, 'different non-zero structure in edge centrality matrix');
end

% make sure betweenness centrality can use an optional edge weight matrix
bc = betweenness_centrality(A,struct('edge_weight',rand(nnz(A),1)));
bc = betweenness_centrality(A);
bc2 = betweenness_centrality(A,struct('edge_weight','matrix'));
if any(bc-bc2)
    error(msgid, 'edge_weight option error');
end
try
    bc = betweenness_centrality(A,struct('edge_weight',rand(2,1)));
    error(msgid, 'betweenness_centrality(weight=rand(2,1)) did not report an error');    
catch
end

%% clustering_coefficients

% Create a clique, where all the clustering coefficients are equal
A = sparse(ones(5));
ccfs = clustering_coefficients(A);
if any(ccfs ~= ccfs(1))
    error(msgid, 'clustering_coefficients failed');
end

%% core_numbers
load('../graphs/kt-7-2.mat');
A = spones(A);
cn = core_numbers(A);
load('../graphs/cores_example.mat');
cn = core_numbers(A);
cn2 = core_numbers(A,struct('unweighted',0));
if any(cn-cn2)
    error(msgid, 'core_numbers failed equivalence test');
end

A = [0 -1 -2; -1 0 -2; -2 -2 0];
cn = core_numbers(sparse(A),struct('unweighted',0));
if any(cn-[-1; -1; -4])
    error(msgid, 'core_numbers failed negative test');
end

%% dominator_tree
load('../graphs/dominator_tree_example.mat');
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= [ 0, 1,  2,  2,  4,   5,  5,  2])
    error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

% graphs from boost example's

A=sparse(13,13);
A(1,2)=1;A(1,3)=1;A(1,4)=1;A(2,5)=1;A(3,2)=1;A(3,5)=1;A(3,6)=1;A(4,7)=1;
A(4,8)=1;A(5,13)=1;A(6,9)=1;A(7,10)=1;A(8,10)=1;A(8,11)=1;A(9,6)=1;
A(9,12)=1;A(10,12)=1;A(11,10)=1;A(12,1)=1;A(12,10)=1;A(13,9)=1;
pred=[0 1 1 1 1 1 4 4 1 1 8 1 5 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

A=sparse(7,7);
A(1,2)=1;A(2,3)=1;A(2,4)=1;A(3,5)=1;A(3,6)=1;A(5,7)=1;A(6,7)=1;A(7,2)=1;
pred=[0 1 2 2 3 3 3 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

A=sparse(13,13);
A(1,2)=1;A(1,3)=1;A(2,4)=1;A(2,7)=1;A(3,5)=1;A(3,8)=1;A(4,6)=1;A(4,7)=1;
A(5,8)=1;A(5,3)=1;A(6,9)=1;A(6,11)=1;A(7,10)=1;A(8,13)=1;A(9,12)=1;
A(10,9)=1;A(11,12)=1;A(12,2)=1;A(12,13)=1;
pred=[0 1 1 2 3 4 2 3 2 7 6 2 1 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

A=sparse(8,8);
A(1,2)=1;A(2,3)=1;A(2,4)=1;A(3,8)=1;A(4,5)=1;A(5,6)=1;A(5,7)=1;A(6,8)=1;
A(7,5)=1;
pred=[0 1 2 2 4 5 5 2 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

A=sparse(8,8);
A(1,2)=1;A(2,3)=1;A(3,4)=1;A(3,5)=1;A(4,3)=1;A(5,6)=1;A(5,7)=1;A(6,8)=1;
A(7,8)=1;
pred=[0 1 2 3 3 5 5 5 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

A=sparse(8,8);
A(1,2)=1;A(1,3)=1;A(2,7)=1;A(3,4)=1;A(3,5)=1;A(4,8)=1;A(6,8)=1;A(7,8)=1;
pred=[0 1 1 3 3 0 2 1 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

A=sparse(14,14);
A(1,2)=1;A(1,14)=1;A(2,3)=1;A(3,4)=1;A(3,8)=1;A(4,5)=1;A(4,6)=1;A(5,7)=1;
A(6,7)=1;A(7,9)=1;A(8,9)=1;A(9,10)=1;A(10,11)=1;A(10,12)=1;A(11,12)=1;
A(12,10)=1;A(12,13)=1;A(13,3)=1;A(13,14)=1;
pred=[0 1 2 3 4 4 4 3 3 9 10 10 12 1 ];
p = lengauer_tarjan_dominator_tree(A,1);
if any(p ~= pred)
   error(msgid, 'lengauer_tarjan_dominator_tree failed test');
end

%% edmonds_maximum_cardinality_matching
load('../graphs/matching_example.mat');
[m,v] = edmonds_maximum_cardinality_matching(A);
if nnz(m)/2 ~= 8
    error(msgid, 'edmonds_maximum_cardinality_matching failed');
end

%% matching
load('../graphs/dfs_example.mat');
try
    [m,v] = matching(A);
    error(msgid,'matching failed');
catch
end
load('../graphs/matching_example.mat');
[m,v] = matching(A);




%% test_dag
% Test the dag_test function, which also tests topological order
A = sparse(6,6);
A(3,6) = 1;
A(1,2) = 1;
A(3,5) = 1;
A(1,4) = 1;
A(2,5) = 1;
A(5,4) = 1;

dag = test_dag(A);
if dag == 0
    error(msgid, 'A dag was not identified as a dag');
end

% Test something that isn't a dag
A = cycle_graph(n);

dag = test_dag(A);
if dag == 1
    error(msgid, 'A cycle was identified as a dag');
end
    
%% test_matching
load('../graphs/matching_example.mat');
[m_max,v] = matching(A);
[m_not_max,v] = matching(A,struct('augmenting_path','none'));
if ~test_matching(A,m_max) || test_matching(A,m_not_max)
    error(msgid, 'test_matching failed');
end    

%% topological_order
% Test the topological order function
n = 10;
A = sparse(1:n-1, 2:n, 1, n, n);
p = topological_order(A);
if any(p - (1:n)')
    error(msgid, 'topological_order failed on simple case');
end