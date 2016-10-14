function caret_metricpermute(P,Oname,varargin)
% function caret_metricpermute(Infiles,Outname,varargin)
% Takes N files with P columns each
% and creates P files with N columns
% INPUTS: 
%   Infiles: character array of input files, prompted if empty
%   Outfile: Name of outfiles. They are written <Outname>_##.img
% _________________________________________________________________
% v.1.0 Joern Diedrichsen 03/06
% jdiedric@bme.jhu.edu

% vararginoptions(varargin,{},{});

if (isempty(P))
    P=spm_get(inf,'*metric',{'Choose metric files ro permute'});
end;
for i=1:length(P)
    M(i)=caret_load(P{i});
    num_cols(i)=M(i).num_cols;
end;
if var(num_cols)>0
    error('Number of columns is not the same');
end;

for file=1:num_cols(i) 
    O.num_cols=size(P,1);
    O.encoding={'BINARY'};
    for col=1:O.num_cols
        O.data(:,col)=M(col).data(:,file);
        O.column_name{col}=[M(col).column_name{file}];
        O.column_color_mapping(col,1:2)=M(col).column_color_mapping(file,1:2);
    end;
    name=sprintf('%s_%2.2d.metric',Oname,file);
    caret_savemetric(name,O);
end;