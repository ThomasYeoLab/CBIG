function caret_savecSPM(filename,cSPM)
% function caret_savecSPM(filename,cSPM)
% saves a statistical cSPM file as a metric file, 
% so that it can be displayed in caret 
% The columns to be saved are: 
% * N columns of the underlying data                   (data_1,...,data_N)
% * P columns of the beta weights                      (B_1,...,B_P)
%     (for the one-sample t-test this would be the mean)                    
% * 1 column of the residual variance                  (ResVar)
% * Number of non-NaN observations per voxel           (N)
% * K statistical values for all contrasts             (Z_1,....Z_K)
% If available: 
% * Normal values corresponding to the p-value of stats (Z_P_1,....,Z_P_K)
% * effect-size values                                (delta_1,...,delta_K)
% ------------------------------------------------------------------
% v.1.0 Joern Diedrichsen 12/04
% jdiedric@bme.jhu.edu
N=size(cSPM.data,2);
P=size(cSPM.b,2);
K=length(cSPM.con);
M.data=[cSPM.data cSPM.b cSPM.ResVar cSPM.N cSPM.con.Z];
where=0;
for n=1:N
    M.column_name{where+n}=sprintf('data_%2.2d',n);
    M.column_color_mapping(where+n,:)=[min(cSPM.data(:)) max(cSPM.data(:))];
end;
where=where+N;
for n=1:P
    M.column_name{where+n}=sprintf('beta_%2.2d',n);
    M.column_color_mapping(where+n,:)=[min(cSPM.b(:,n)) max(cSPM.b(n,:))];
end;
where=where+P;
M.column_name{where+1}=sprintf('ResVar',n);
M.column_color_mapping(where+1,:)=[min(cSPM.ResVar) max(cSPM.ResVar)];
M.column_name{where+2}=sprintf('N',n);
M.column_color_mapping(where+2,:)=[min(cSPM.N) max(cSPM.N)];
where=where+2;
for n=1:K
    M.column_name{where+n}=sprintf('Z_%d',n);
    M.column_color_mapping(where+n,:)=[min(cSPM.con(n).Z) max(cSPM.con(n).Z)];
end;
where=where+K;
if (isfield(cSPM.con(1),'Z_P'))
    for n=1:K
        M.column_name{where+n}=sprintf('Z_P_%d',n);
        M.column_color_mapping(where+n,:)=[min(cSPM.con(n).Z_P) max(cSPM.con(n).Z_P)];
    end;
    M.data=[M.data cSPM.con(n).Z_P];
    where=where+K;
end;
if (isfield(cSPM.con(1),'delta'))
    for n=1:K
        M.column_name{where+n}=sprintf('delta_%d',n);
        M.column_color_mapping(where+n,:)=[min(cSPM.con(n).delta) max(cSPM.con(n).delta)];
    end;
    M.data=[M.data cSPM.con(n).delta];
end;
M.encoding={'BINARY'};
caret_savemetric(filename,M);
