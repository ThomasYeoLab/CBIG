function test_layouts

%%
msgid = 'matlab_bgl:test_layouts';

%% circle_graph_layout
X = circle_graph_layout(sparse(4,4));
Y = [1 0; 0 1; -1 0; 0 -1];
if norm(X-Y)>5e-6, error(msgid, 'circle_graph_layout(4) bad coords'); 
end

%% fruchterman_reingold_force_directed_layout
G = cycle_graph(10,struct('directed',0));
X = fruchterman_reingold_force_directed_layout(G);
G = cycle_graph(10);
try
  G = fruchterman_reingold_force_directed_layout(G);
  error(msgid, 'fruchterman_reingold failed to throw on a directed graph');
catch end
G = grid_graph(5,6);
X = fruchterman_reingold_force_directed_layout(G);
for i=0:10, X = fruchterman_reingold_force_directed_layout(sparse(i,i)); end
G = cycle_graph(10,struct('directed',0));
X1 = fruchterman_reingold_force_directed_layout(G,...
  'progressive',circle_graph_layout(G));
X2 = fruchterman_reingold_force_directed_layout(G,...
  struct('progressive',circle_graph_layout(G)));
if norm(X1-X2,'fro')>1e-10, error(msgid, 'fruchterman_reingold options failed'); end



%% gursoy_atun_layout
check_layout = @(X) all(all(isfinite(X)));
for i=0:10, X=gursoy_atun_layout(sparse(i,i)); assert(check_layout(X)); end
for i=0:10, X=gursoy_atun_layout(grid_graph(i,i)); assert(check_layout(X)); end
for i=0:10, for j=2:10,
        X=gursoy_atun_layout(wheel_graph(i),'topology',sprintf('cube%i',j));
        assert(check_layout(X)); 
        X=gursoy_atun_layout(wheel_graph(i),'topology',sprintf('ball%i',j));
        assert(check_layout(X));         
end, end

try
    X=gursoy_atun_layout(sparse(10,10),'topology','ball');
    error(msgid,'gursoy_atun_layout did not throw on invalid topology(ball)');
catch end

try
    X=gursoy_atun_layout(sparse(10,10),'topology','cube');
    error(msgid,'gursoy_atun_layout did not throw on invalid topology(cube)');
catch end

%% kamada_kawai_spring_layout
% TODO: Expand these test cases
try
  X = kamada_kawai_spring_layout(sparse(2,2));
  error(msgid,'kamada_kawai_spring_layout didn''t fail on empty graph');
catch
end
X = kamada_kawai_spring_layout(grid_graph(1,2));
for i=0:10, X=kamada_kawai_spring_layout(grid_graph(i,i)); assert(check_layout(X)); end


%% random_graph_layout
rand('state',0);
X = random_graph_layout(sparse(4,4));
if ~isequal(size(X),[4,2])
  error(msgid,'random_graph_layout(4) wrong output size');
end
if any(any(X<0)) || any(any(X>1))
  error(msgid,'random_graph_layout(4) wrong output area');
end
X = random_graph_layout(sparse(1500,1500),[-1e5,-1e5,1e5,1e5]);
if ~isequal(size(X),[1500,2])
  error(msgid,'random_graph_layout(1500) wrong output size');
end
if any(any(X<-1e5)) || any(any(X>1e5))
  error(msgid,'random_graph_layout(1500) wrong output area');
end

X = random_graph_layout(sparse(256,256),int32([0 0 0 0]));
if ~all(all(X==0))
  error(msgid,'random_graph_layout(256) invalid output');
end
