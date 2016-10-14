function caret_metriccalc(P,Oname,form,varargin)
% function MAP=caret_metriccalc(Infiles,Out,form,varargin)
% Calculator for metric files. 
% All metric files have to have the same number of columns 
% INPUTS:
%    Infiles: input metric files 
%    Out: name of output files 
%    form: string to be evaluated: 
%         'i#' is metric file #n
%      in matrix mode: 
%      all metric files are put into a matrix called X 
%    varargin: 'matrix'
% EXAMPLE:
% caret_metriccalc({'files1.metric','file2.metric'},'diff.metric','i1-i2');
%               Calculates difference between file1 and file2, column by column
% caret_metriccalc({'files1.metric','file2.metric',...},'sum.metric','sum(X)','matrix');
% _________________________________________________________________
% v.1.0 Joern Diedrichsen 01/05
% jdiedric@bme.jhu.edu
c=1;
is_matrix=0;
if (iscell(P))
    P=char(P);
end;
while c<=length(varargin)
    switch (varargin{c})
        case 'matrix'
            is_matrix=1;   
            c=c+1;
        otherwise 
            error(['Unknown option:' varargin{c}]);
    end;
end;
if (isempty(P))
    P=spm_get(inf,'*metric',{'Choose metric files to do calculation on'});
end;
for i=1:size(P,1)
    M(i)=caret_load(deblank(P(i,:)));
    num_cols(i)=M(i).num_cols;
end;
if var(num_cols)>0
    error('Number of columns is not the same');
end;

O.num_cols=num_cols(i);
O.encoding={'BINARY'};
for c=1:num_cols(i)
    for i=1:size(P,1)
        if (is_matrix)
            X(i,:)=M(i).data(:,c)';
        else
            eval(sprintf('i%d = transpose(M(i).data(:,c));',i));
        end;
    
    end;
    O.data(:,c)=eval(form)';
    O.column_name{c}=[form ' on ' M(1).column_name{c}];
    O.column_color_mapping(c,1:2)=[-3 3];
end;
if (isempty(O))
    [Oname,odir]=uiputfile('*.metric','Save metric as');
end;
caret_savemetric(Oname,O);
