function cSPM=caret_getcSPM(type,varargin)
% function MAP=caret_getcSPM(type,varargin)
% generates a statistical map from columns of metric files 
%   type: Type of linear model 
%       'onesample_t':  tests the Hypothesis that the mean effect is bigger
%                       than zero
%                       the only implemented type of map right now is the
%       'regression':   simple or multiple regression  
%                       'regression',X 
%                       'regression',X,'no_intercept'
%                       where X is a regression matrix 
%                       The intercept term will be added, unless otherwise
%                       specified
%   varargin: 
%       'data',filename,columns: 
%                       Specifies where the data comes from. To concatinate
%                       data from different files, just repeat the argument
%       'Z_P':          adds to each contrast a Normal value, corresponding to the
%                       p-value of the contrast
%       'delta':        add to each contrast the effect-size (mean/SD) for
%                       comparision of effects with different number of
%                       subjects 
%        'maskthreshold',n: only computes the statistics when N>=n subject are 
%                       measured at this node 
% _________________________________________________________________
% v.1.1 Joern Diedrichsen 01/04
% jdiedric@bme.jhu.edu
if(nargin<1 | isempty(type)) 
    type='onesample_t'
end; 
cSPM.data=[];
cSPM.N=0;
switch (type)
    case 'onesample_t'
        c=1;
    case 'regression'
        X=varargin{1};
        c=2;
    otherwise 
        error(['Unknown type: ' type]);
end;
maskthreshold=0;
z_p=0;
delta=0;
is_intercept=1;
while c<=length(varargin)
    switch (varargin{c})
        case 'no_intercept'
            is_intercept=0;    
        case 'data'
            try 
                filename=varargin{c+1};
                columns=varargin{c+2};
            catch 
                error('Syntax: caret_getcSPM(type,"data",filename,columns,...)');    
            end;
            M=caret_load(filename);
            if (isempty(columns))
                columns=[1:M.num_cols];
            end;
            cSPM.data=[cSPM.data M.data(:,columns)];
            c=c+3;
        case 'Z_P'             % Compute and write out the z-value corresponding to the P-value
            z_p=1;
            c=c+1;
        case 'delta'            % Compute and write out effect-size
            delta=1;
            c=c+1;
        case 'maskthreshold'
            maskthreshold=varargin{c+1};
            c=c+2;
        otherwise 
            error(['Unknown option:' varargin{c}]);
    end;
end;
N=size(cSPM.data,2);
j=find(cSPM.data==0.0);
cSPM.data(j)=NaN;
cSPM.N=sum(~isnan(cSPM.data'))';
switch (type)
    case 'onesample_t'
        X=ones(N,1);cSPM.X=X;
        cSPM.b=nanmean(M.data')';
        cSPM.ResVar=(nanstd(cSPM.data')').^2;
        cSPM.ResVar(cSPM.ResVar==0)=NaN;
        cSPM.bvar=inv(X'*X);
        cSPM.P=1;
        cSPM.title='one-sample t-test: mean>0';
    case 'regression'
        if (is_intercept)
            X=[ones(N,1) X];
        end;
        cSPM.X=X;
        iX=inv(X'*X);
        cSPM.b=(iX*X'*cSPM.data')';
        res=(cSPM.data'-cSPM.X*cSPM.b')';
        cSPM.bvar=iX;
        cSPM.P=rank(X);
        cSPM.ResVar=sum(res.^2,2)./(cSPM.N-cSPM.P);
        cSPM.title='regression t-test: b>0';
end;

for i=1:size(X,2)
    cSPM.con(i).df=[1 max(cSPM.N)-cSPM.P];
    cSPM.con(i).STAT='T';
    c=zeros(cSPM.P,1);c(i)=1;
    cSPM.con(i).Z=(c'*cSPM.b'./sqrt((c'*cSPM.bvar*c)*cSPM.ResVar'))';  
    if (maskthreshold>0)
        cSPM.con(i).Z(cSPM.N<maskthreshold)=NaN;
    end;
    cSPM.con(i).Z(isnan(cSPM.con(i).Z))=0;
    if (z_p)
        cSPM.con(i).Z_P=norminv(tcdf(cSPM.con(i).Z,cSPM.con(i).df(2)));
    end;
    if (delta)
        cSPM.con(i).delta=c'*cSPM.b(:,i)./sqrt(cSPM.ResVar);
    end;
end;

