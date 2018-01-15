function test_main
%% Implement a test suite for matlab_bgl
% Call this function to test the MatlabBGL package.

%% Setup some parameters

msgid = 'matlab_bgl:test_main';

%%
fprintf('testing trivial input...\n');
test_trivial

%%
fprintf('testing code examples...\n');
test_examples

%% 
fprintf('testing old test script...\n');
test_all

%%
fprintf('testing for regressions...\n');
rtest_all

%%
fprintf('testing component algs...\n');
test_components

%% 
fprintf('testing layout algs...\n');
test_layouts

%% 
fprintf('testing planar graph algs...\n');
test_planar

%%
fprintf('testing search algs...\n');
test_searches

%%
fprintf('testing shorest path algs...\n');
test_shortest_paths

%% 
fprintf('testing spanning tree algs...\n');
test_spanning_trees

%% 
fprintf('testing statistics algs..\n');
test_statistics

%%
fprintf('testing utility funcs...\n');

%% cycle_graph
G = cycle_graph(10,struct('directed',0));
G = cycle_graph(10,'directed',0);
[A,xy] = cycle_graph(0); assert(all(size(A)==[0 0]));
[A,xy] = cycle_graph(1); assert(all(size(A)==[1 1])); assert(nnz(A)==1);
for i=0:10, [A,xy] = cycle_graph(i); end
for i=50:50:250, [A,xy] = cycle_graph(i); end

%% edge_weight_vector
n = 8; u = 1; v = 2;
E = [1:n 2:n 1; 2:n 1 1:n]';
w = [1 zeros(1,n-1) 1 zeros(1,n-1)]';
A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
[d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));
if d(v) ~= 0
    error(msgid, 'edge_weight_vector failed');
end

% remove the edge between node 1 and edge 8 to test non-symmetric in As
A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
As(1,8) = 0;
[d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));
if d(v) ~= 1 || any(d>1) 
    error(msgid, 'edge_weight_vector failed');
end

% make the weights non-symmetric to test non-symmetry in A
A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
A(1,2) = 2; 
A(1,8) = 3;
A(2,3) = 4;
[d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));
if d(v) ~= 2 || any(d(d>2)~=3)
    error(msgid, 'edge_weight_vector failed');
end

% test non-symmetric A and non-symmetric As
A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
As(1,8) = 0;
A(1,2) = 2; 
A(2,3) = 4;
[d pred] = shortest_paths(As,u,struct('edge_weight',edge_weight_vector(As,A)));
if d(v) ~= 2 || any(d(d>2)~=6)
    error(msgid, 'edge_weight_vector failed');
end

% make sure it works with pre-transposed matrices, repeat all the test
% cases

A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
[d pred] = shortest_paths(As',u,...
    struct('edge_weight',edge_weight_vector(As',A',struct('istrans',1)),...
        'istrans',1));
if d(v) ~= 0
    error(msgid, 'edge_weight_vector failed');
end

A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
As(1,8) = 0;
[d pred] = shortest_paths(As',u,...
    struct('edge_weight',edge_weight_vector(As',A',struct('istrans',1)),...
        'istrans',1));
if d(v) ~= 1
    error(msgid, 'edge_weight_vector failed');
end

A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
A(1,2) = 2; 
A(1,8) = 3;
A(2,3) = 4;
[d pred] = shortest_paths(As',u,...
    struct('edge_weight',edge_weight_vector(As',A',struct('istrans',1)),...
        'istrans',1));
if d(v) ~= 2 || any(d(d>2)~=3)
    error(msgid, 'edge_weight_vector failed');
end

A = sparse(E(:,1), E(:,2), w, n, n); % create weighted sparse matrix
As = sparse(E(:,1), E(:,2), true, n, n); % create structural sparse matrix
As(1,8) = 0;
A(1,2) = 2; 
A(2,3) = 4;
[d pred] = shortest_paths(As',u,...
    struct('edge_weight',edge_weight_vector(As',A',struct('istrans',1)),...
        'istrans',1));
if d(v) ~= 2 || any(d(d>2)~=6)
    error(msgid, 'edge_weight_vector failed');
end



%% max_flow

%% pred_from_path

% Create a line graph
n = 10;
A = sparse(1:n-1,2:n,1,n,n);
A = A+A';

% Compute BFS and test pred_from_path
[d dt pred] = bfs(A,1);
path = path_from_pred(pred,n);
if any(path ~= 1:n)
    error(msgid, 'path_from_pred failed');
end




%% star_graph
[A,xy] = star_graph(0); assert(all(size(A)==[0 0]));
[A,xy] = star_graph(1); assert(all(size(A)==[1 1])); assert(nnz(A)==0);
for i=0:10, [A,xy] = star_graph(i); end
for i=50:50:250, [A,xy] = star_graph(i); end

%% wheel_graph
[A,xy] = wheel_graph(0); assert(all(size(A)==[0 0]));
[A,xy] = wheel_graph(1); assert(all(size(A)==[1 1])); assert(nnz(A)==0);
for i=0:10, [A,xy] = wheel_graph(i); end
for i=50:50:250, [A,xy] = wheel_graph(i); end

%%
% ***** end test_main *****
end

