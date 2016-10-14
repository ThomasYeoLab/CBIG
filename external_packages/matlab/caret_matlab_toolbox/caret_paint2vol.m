function T=caret_paint2vol(paint_name,coord_name,target_name,varargin);
%  function T=caret_paint2vol(roi_name,image_names,vararin);
%       Translates a paint file into a volume file  
% -------------------------------------------------------------------------
% v.1.0 Joern Diedrichsen jdiedric@bme.jhu.edu
c=1;
if (nargin<1 | isempty(paint_name)) 
    paint_name=spm_get(1,'.paint','Select paint file');
end;
if (nargin<2 | isempty(coord_name)) 
    coord_name=spm_get(1,'.coord','Select Fiducial file');
end;
if (nargin<3 | isempty(target_name)) 
    target_name=spm_get(1,'.img','Space defining image');
end;

Width=1;
Convolution_kernel=ones([Width Width Width]);
PAINT=caret_load(paint_name);
COORD=caret_load(coord_name);
V=spm_vol(target_name); 
for s=1:PAINT.num_cols 
    PaintVol=V;
    PaintVol.fname=[PAINT.column_name{s} '.img'];
    X=zeros(V.dim(1:3));
    indx=find(PAINT.data(:,s)>0);
    PIXEL=[COORD.data(indx,:) ones(length(indx),1)];
    PIXEL_tr=round(inv(V.mat)*PIXEL')';
    for i=1:length(indx)
        X(PIXEL_tr(i,1),PIXEL_tr(i,2),PIXEL_tr(i,3))=PAINT.data(indx(i),s);
    end;
    X=convn(X,Convolution_kernel,'same');
    X(X>0)=1;
    PaintVol=spm_create_vol(PaintVol);
    spm_write_vol(PaintVol,X);
end;
