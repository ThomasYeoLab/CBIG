function index=caret_clusters(Surface,i);
% Clusters a surface/Subsurface 
global gS;
global visited;
gS=Surface;
cluster=1;
if (nargin<2)
    visited=zeros(gS.num_nodes,1);
else 
    visited=zeros(gS.num_nodes,1);
    visited(find(i==0))=-1;
end;

isfinished=0;
while (~isfinished)
    N=find(visited==0);
    if (length(N)==0)
        isfinished=1;
    else
        fprintf('.');
        caret_bounce(N(1),cluster);
        cluster=cluster+1;
    end;
end;
index=visited;
index(index==-1)=0;
fprintf('\n');
    

function caret_bounce(n,cluster);
global gS;
global visited;
visited(n)=cluster;
for N=gS.Nodes.Neighbor{n}
    if (visited(N)==0)
        caret_bounce(N,cluster);    
    end;
end;
