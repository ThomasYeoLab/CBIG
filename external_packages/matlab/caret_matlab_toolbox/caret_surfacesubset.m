function SubS=caret_surfacesubset(S,i);
% ---------------------------------------------
% function SubS=caret_surfacesubset(S,i);
% Transfer Nodes, generate the new node numbers 
add=find(i>0);
new_nodenumber=zeros(S.num_nodes,1);
new_nodenumber(add)=[1:length(add)]';
SubS.num_nodes=length(add);
% Transfer the Nodes
SubS.Nodes.data=S.Nodes.data(add,:);
SubS.Nodes.area=S.Nodes.area(add,:);

% ---------------------------------------------
% Transfer Tiles, update the node numbers 
T=zeros(S.num_tiles,3);
for j=1:3 
    T(:,j)=i(S.Tiles.data(:,j));
end;
addT=find(all(T')');
SubS.Tiles.data=S.Tiles.data(addT,:);
for j=1:3 
    SubS.Tiles.data(:,j)=new_nodenumber(S.Tiles.data(addT,j));
end;

if(isfield(S.Tiles,'area'))
    SubS.Tiles.area=S.Tiles.area(addT,:);
end;
if(isfield(S.Tiles,'edgelength'))
    SubS.Tiles.edgelength=S.Tiles.edgelength(addT,:);
end;
SubS.num_tiles=length(addT);
SubS.sub_index=add;
