function []=CBIG_gwMRF_write_cifti_from_annot(lh_annot,rh_annot,output_name,nparcels,lh_fslr_orig,rh_fslr_orig)

% [] = CBIG_gwMRF_write_cifti_from_annot(base_folder)
%
%This script writes files from freesurfer annot format into HCP cifti
%format. It has dependencies on CBIG, Workbench and Fieldtrip. Further it
%also needs to load in an 'common_cifti.mat' file.

%input 
% -lh_annot = annot file containing the orginal parcellation in freesurfer
% format
% -rh_annot = annot file containing the orginal parcellation in freesurfer
% format
% -output_name = Naming of final cifti file 
% -nparcels =  Number of parcels in each of the hemispheres
% -lh_fslr_orig = left hemisphere vector of same parcellation which has
% already been projected to fslr32k
% -rh_fslr_orig = right hemisphere vector of same parcellation which has
% already been projected to fslr32k
%
%output is generated as files
%
%example
%CBIG_gwMRF_write_cifti_from_annot(lh_annot,rh_annot,output_name,500,lh_fslr32k,rh_fslr32k)

% Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
load(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'Code', 'input', 'common_cifti.mat'));
[a,lh,lh_s]=read_annotation(lh_annot);
[a,rh,rh_s]=read_annotation(rh_annot);

lh_new=zeros(size(lh));
rh_new=zeros(size(rh));
for i=1:nparcels
    lh_new(lh==lh_s.table(i+1,5))=i;
    rh_new(rh==rh_s.table(i+1,5))=i;
end

fileID=fopen([output_name,'_info.txt'],'w');

for i=2:nparcels+1
        fprintf(fileID,'%s\n',lh_s.struct_names{i});
        fprintf(fileID,'%i %i %i %i %i\n',i-1,lh_s.table(i,1),lh_s.table(i,2),lh_s.table(i,3),255);
end
for i=2:nparcels+1
        fprintf(fileID,'%s\n',rh_s.struct_names{i});
        fprintf(fileID,'%i %i %i %i %i\n',nparcels+i-1,rh_s.table(i,1),rh_s.table(i,2),rh_s.table(i,3),255);
end
fclose(fileID);

system('rm -rf ~/temp123/')
[lh_fslr,rh_fslr]=CBIG_project_fsaverage2fsLR(lh_new,rh_new,'fsaverage6','label','~/temp123/','20160827');
% Exclude medial wall vertices defined by HCP
lh_mesh_fslr_32k=CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','medialwall.annot');
rh_mesh_fslr_32k=CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','medialwall.annot');
lh_fslr(lh_mesh_fslr_32k.MARS_label == 1) = 0;
rh_fslr(rh_mesh_fslr_32k.MARS_label == 1) = 0;

%% to avoid that interpolation errors cause differences 
lh_label_max_match=zeros(size(lh_fslr));
rh_label_max_match=zeros(size(rh_fslr));
for i=1:max(lh_fslr_orig);
     lh_label_max_match(lh_fslr_orig==i)=mode(lh_fslr(lh_fslr_orig==i));
end
for i=1:max(rh_fslr_orig);
  rh_label_max_match(rh_fslr_orig==i)=mode(rh_fslr(rh_fslr_orig==i));
end
rh_label_max_match=rh_label_max_match+nparcels;

%%% remove isolated components
lh_mesh_fslr_32k=CBIG_read_fslr_surface('lh','fs_LR_32k','inflated','aparc.annot');
lh_label_max_match = CBIG_RemoveIsolatedSurfaceComponents(lh_mesh_fslr_32k,lh_label_max_match, 3);
rh_mesh_fslr_32k=CBIG_read_fslr_surface('rh','fs_LR_32k','inflated','aparc.annot');
rh_label_max_match = CBIG_RemoveIsolatedSurfaceComponents(rh_mesh_fslr_32k,rh_label_max_match, 3);
rh_label_max_match(rh_label_max_match==min(rh_label_max_match))=0;
fslr.parcels=double([lh_label_max_match',rh_label_max_match'])';
fslr.parcelslabel=[lh_s.struct_names(2:nparcels+1);rh_s.struct_names(2:nparcels+1)];
fslr.pos=pos;
fslr.brainstructure=brainstructure;
fslr.brainstructurelabel=brainstructurelabel';
fslr.dimord='pos';
ft_write_cifti(output_name,fslr,'parcellation','parcels','parameter','parcels');
system(['wb_command -cifti-label-import ',output_name,'.dscalar.nii ',output_name,'_info.txt ',...
    output_name,'.dlabel.nii'])
