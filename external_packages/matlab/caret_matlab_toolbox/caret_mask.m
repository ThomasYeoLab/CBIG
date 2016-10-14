function caret_mask(P,mask,column,expression,varargin)
% function caret_mask(P,M,column,expression)
% Does the masking of a metric file (setting out of mask to zero)
% mask is the expression with x for the column
% EXAMPLE:
% caret_mask('myfile.metric','mymask.metric',5,x>5);
%   Sets all values in the metric file 'myfile' to zero, for which 
%   the 5th column of mymask.metric is smaller 5.
% _________________________________________________________________
% v.1.0 Joern Diedrichsen 01/05
% jdiedric@jhu.edu
c=1;
while c<=length(varargin)
    switch (varargin{c})
        otherwise 
            error(['Unknown option:' varargin{c}]);
    end;
end;
if (isempty(P))
    P=spm_get(inf,'*metric',{'Choose metric to mask'});
end;
if (~iscell(P))
    P={P};
end;
for i=1:length(P)
    M(i)=caret_load(P{i});
    num_cols(i)=M(i).num_cols;
end;

if (isempty(mask))
    mask=spm_get(1,'*metric',{'Choose mask-file'});
end;
if (iscell(mask))
    mask=mask{1};
end;
MASK=caret_load(mask);
x=MASK.data(:,column);
indx=eval(expression);

for i=1:length(P)
    for j=1:M(i).num_cols;
        M(i).data(indx==0,j)=0;
    end;
    caret_savemetric(P{i},M(i));
end;

