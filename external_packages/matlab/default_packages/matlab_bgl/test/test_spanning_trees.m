function test_spanning_trees

msgid= 'matlab_bgl:test_spanning_trees';

%% mst
% make sure mst('kruskal','root') fails
load('../graphs/clr-24-1.mat');

try
    T = mst(A,struct('algname','kruskal','root',5));
    error(msgid,'mst(kruskal,root) did not report an error');
catch
end

%% prim_mst
load('../graphs/clr-24-1.mat');
% change the graph to make the result unique
A(2,3) = 9; A(3,2) = 9;
T = prim_mst(A); % root the tree at vertex 
Ttrueijv = [
     2     8     1     4     6     9     3     5     4     3     7     6     8     1     7     3
     1     1     2     3     3     3     4     4     5     6     6     7     7     8     8     9
     4     8     4     7     4     2     7     9     9     4     2     2     1     8     1     2 ];
Ttrue = sparse(Ttrueijv(1,:),  Ttrueijv(2,:), Ttrueijv(3,:), 9,9);
if nnz(T - Ttrue) ~= 0
    error(msgid, 'prim_mst failed');
end

% unfortunately, I can't deterministically check if the 
% rooted option works because which tree the algorithm
% will pick is implementation dependent and not guaranteed
load('../graphs/clr-24-1.mat');
% this example should work, however...
T1 = prim_mst(A); % root the tree at vertex 
T2 = prim_mst(A,struct('root',5));
T1T2diff = sparse(9,9); T1T2diff(2,3) = -8; T1T2diff(1,8) = 8;
if nnz(triu(T1-T2)-T1T2diff) ~= 0
    error(msgid, 'prim_mst failed rooted test');
end
