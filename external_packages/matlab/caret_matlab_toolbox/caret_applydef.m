function M=caret_applydef(infilename,deformation_map,outfilename)
% function M=caret_applydef(infilename,outfilename,deformation_map)
% Applies a spm-deformation map to a caret coord file
% 
%---------------------------------------------------
% v.1.0 Joern Diedrichsen 2/06 
% jdiedric@bme.jhu.edu
if(nargin<1)
    [file,path]=uigetfile('*.coord','coord file to deform');
    infilename=[path file];
end;

if(nargin<2)
    [file,path]=uigetfile('*.img','deformation map');
    deformation_map=[path file];
end;

COORD=caret_load(infilename); 
for i=1:3 
    Vdef(i)=spm_vol(sprintf('%s,%d',deformation_map,i));
end;
[x,y,z]=spmj_affine_transform(COORD.data(:,1),COORD.data(:,2),COORD.data(:,3),inv(Vdef(1).mat));
Nx=spm_sample_vol(Vdef(1),x,y,z,1);
Ny=spm_sample_vol(Vdef(2),x,y,z,1);
Nz=spm_sample_vol(Vdef(3),x,y,z,1);
COORD.data=[Nx Ny Nz];
if(nargin<3)
    [file,path]=uiputfile('*.coord','Outfile');
    outfilename=[path file];
end;

caret_save(outfilename,COORD);
