function T=caret_regionstat(roi_name,image_names,varargin);
%  function T=caret_regionstat(roi_names,image_names,vararin);
%       Extracts the values of the images over the regions 
%       defined in the columns of paint files 
% -------------------------------------------------------------------------
% v.1.0 Joern Diedrichsen jdiedric@bme.jhu.edu
c=1;
T=[];T.region=[];T.region_name={};T.image=[];T.image_name={};T.column=[];T.column_name={};
T.reg_size=[];T.reg_mean=[];

if (nargin<1 | isempty(roi_name)) 
    roi_name=spm_get(1,'.paint','Select regions');
end;
if (nargin<2 | isempty(image_names))
    image_names=spm_get([1:1000],'.metric','Select files');
end;

R=caret_load(roi_name);
for i=1:size(image_names,1)
    M=caret_load(deblank(image_names(i,:)));
    fprintf('Metric: %s\n',deblank(image_names(i,:)));
    for r=1:R.num_cols
        indx=find(R.data(:,r)>0);
        for c=1:M.num_cols
            T.region(end+1,1)=r;
            T.region_name{end+1,1}=R.column_name{r};
            T.image(end+1,1)=i;
            T.image_name{end+1,1}=image_names(i,:);
            T.column(end+1,1)=c;
            T.column_name{end+1,1}=M.column_name{c};
            T.reg_size(end+1,1)=length(indx);
            d=M.data(indx,c); 
            if sum(~isnan(d))==0
                T.reg_mean(end+1,1)=NaN;
            else 
                T.reg_mean(end+1,1)=nanmean(d);
            end;
        end;    
    end;
end;
