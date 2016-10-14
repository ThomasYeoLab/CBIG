function data=caret_threshold(data,S,u,k,columns);
% function data=caret_threshold(data,S,u,k,columns);
% thresholds the map M on surface S at a height-threshold u and area
% threshold k (in mm2). All other values are set to 0.
% data: N*p data matrix to be thresholded (M.data or cSPM.con.Z) 
% S: surface structure: Neighborhood has to be calculated 
% u: Height threshold. if u is a scalar, it selects y>u 
%       if u is a 1*2 vector, then it select y<u(1) or y>u(2) 
%       for only a negative threshold, set u=[-x inf];
% k: cluster size threshold in mm^2 
% ---------------------------------------------------------------
% v.1.0 Joern Diedrichsen 05/01/03 
[N,p]=size(data);
if (nargin<5) 
    columns=[1:p];
end;
if (N~=S.num_nodes) 
    error('Data structure has different numbers of nodes than surface');
end;
if (length(u)==1)
    indx=data>u
elseif (length(u)==2)
    indx=data<u(1) | data>u(2);
else 
    error('u has to be either 1 or 2 length');
end;

% Kill data below threshold
for c=columns
    data(indx(:,c)==0,c)=0;
end;

% Kill clusters below threshold 
if (k>0) 
    for c=columns
        A=caret_clusters(S,indx(:,c));
        a=unique(A);
        for j=1:length(a)-1
            i=find(A==j);
            area=sum(S.Nodes.area(i));
            K(j,1)=area;
            N(j,1)=length(i);
            if (area<k) 
                data(i,c)=0;   
            end;
        end;
    end;
end;
