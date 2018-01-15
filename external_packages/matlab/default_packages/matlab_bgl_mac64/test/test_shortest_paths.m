function test_shortest_paths

%%
msgid = 'matlab_bgl:test_shortest_paths';

%% all_shortest_paths
Dtrue = [0 1 -3 2 -4; 3 0 -4 1 -1; 7 4 0 5 3; 2 -1 -5 0 -2; 8 5 1 6 0];
load('../graphs/clr-26-1.mat');
D = all_shortest_paths(A);
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths returned an incorrect distance matrix');
end
D = all_shortest_paths(A,struct('algname','johnson'));
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths(johnson) returned an incorrect distance matrix');
end
D = all_shortest_paths(A,struct('algname','floyd_warshall'));
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths(floyd_warshall) returned an incorrect distance matrix');
end
At = A';
D = all_shortest_paths(At,struct('istrans',1));
if any(any(D' ~= Dtrue))
    error(msgid, 'all_shortest_paths(istrans=1) returned an incorrect distance matrix');
end
% test non-reachable vertex
A = sparse([1 1; 0 1]);
Dtrue = [0 1; 5 0];
D = all_shortest_paths(A,struct('inf',5,'algname','johnson'));
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths(johnson,inf=5) returned an incorrect distance matrix');
end
D = all_shortest_paths(A,struct('inf',5,'algname','floyd_warshall'));
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths(floyd_warshall,inf=5) returned an incorrect distance matrix');
end
% test edge weighted graph
Dtrue = [0 1 -3 2 -4; 3 0 -4 1 -1; 7 4 0 5 3; 2 -1 -5 0 -2; 8 5 1 6 0];
load('../graphs/clr-26-1.mat');
D = all_shortest_paths(A,struct('edge_weight','matrix'));
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths(weight=matrix) returned an incorrect distance matrix');
end
v=nonzeros(A');
D = all_shortest_paths(spones(A),struct('edge_weight',v));
if any(any(D ~= Dtrue))
    error(msgid, 'all_shortest_paths(weight=matrix) returned an incorrect distance matrix');
end
try
    bc = all_shortest_paths(A,struct('edge_weight',rand(2,1)));
    error(msgid, 'all_shortest_paths(weight=rand(2,1)) did not report an error');    
catch
end

% test predecessor matrix
[D P] = all_shortest_paths(A,struct('algname','floyd_warshall'));
for i=1:size(A,1)
    [d p] = shortest_paths(A,i);
    if any(D(i,:)~=d'), error(msgid,'all_shortest_paths(floyd_warshall) returned incorrect distance'); end
    if any(P(i,:)~=p), error(msgid,'all_shortest_paths(floyd_warshall) returned incorrect predecessor'); end    
    % the following command should always work.
    for j=1:size(A,1)
        p=[]; while j~=0, p(end+1)=j; j=P(i,j); end; p=fliplr(p);
    end
end

%% shortest_paths