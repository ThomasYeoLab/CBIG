function caret_avrgPALSsurface(Side);
dir='F:/Atlas_templates/Mapping_files';
if Side=='L'
    topo_name='Human.sphere_6.LEFT_HEM.73730.topo';
elseif Side=='R'
    topo_name='Human.sphere_6.RIGHT_HEM.73730.topo';
end;    

for i=1:12
    if (i<7)
        Gender='F';
    else
        Gender='M';
    end;
    coord_name=sprintf('%s/Human_Buck_Case%d.%s.%s.RegToPALS_B12.LR.FIDUCIAL_SPM2.73730.coord',dir,i,Side,Gender);
    S{i}=caret_getsurface(coord_name,topo_name);
    S{i}=caret_calcarea(S{i});
end;
for i=1:12
    MEAN=S{2};
    Nodes.data(:,:,i)=S{i}.Nodes.data;
    Nodes.ath=mean(Tiles.edgelength,3);
    MEAN.Tiles.area=mean(Tiles.area,2);
    MEAN.Edges.rea(:,i)=S{i}.Nodes.area;
    Tiles.edgelength(:,:,i)=S{i}.Tiles.edgelength;
    Tiles.area(:,i)=S{i}.Tiles.area;
    Edges.Length(:,i)=S{i}.Edges.Length;
end;
MEAN.Nodes.data=mean(Nodes.data,3);
MEAN.Nodes.area=mean(Nodes.area,2);
MEAN.Tiles.endgelengLength=mean(Edges.Length,2);
MEAN.A(3)=sum(MEAN.Nodes.area);
M=caret_struct('metric','data',[Nodes.area MEAN.Nodes.area]);
M.column_name{end}='mean';
caret_save(sprintf('PALS.Avrg.%s.SPM2.Area.metric',Side),M);
S=MEAN;
save(sprintf('PALS.Avrg.%s.SPM2.mat',Side),'S');