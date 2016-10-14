function S=caret_getsurface(coord_file,topo_file);
% function S=caret_getsurface(coord_file,topo_file);
%   Loads a coordinate File and a topology file and 
%   Creates a surface structure, for later usage
%       S.coord: holds coordinate file name and header
%       S.topo: hols topology file name and header
%       S.num_nodes: number of Nodes
%       S.num_tiles: Number of Tiles 
%       S.Nodes: Data for Nodes
%       S.Tiles: Data for Tiles 
%---------------------------------------------------
% v.1.1 Joern Diedrichsen 02/05 
% jdiedric@bme.jhu.edu

if (nargin<1)
    [file,path]=uigetfile('*.coord','get coordinate file');
    S.coord.name=[path file];
else
    S.coord.name=coord_file;
end;
X=caret_load(S.coord.name);
if (nargin<2)
    [file,path]=uigetfile('*.topo','get topology file');
    S.topo.name=[path file];
else 
    S.topo.name=topo_file;
end;
T=caret_load(S.topo.name);
S.num_nodes=X.num_nodes;
S.num_tiles=T.num_tiles;
S.Nodes.data=X.data;
S.Tiles.data=T.data;
if (isfield(T,'Neighbor'))
    S.Nodes.Neighbor=T.Neighbor;
end;
S.topo.header=T.header;
S.coord.header=X.header;