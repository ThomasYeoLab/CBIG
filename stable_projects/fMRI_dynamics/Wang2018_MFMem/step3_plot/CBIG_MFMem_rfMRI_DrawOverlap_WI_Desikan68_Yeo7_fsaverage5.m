function CBIG_MFMem_rfMRI_DrawOverlap_WI_Desikan68_Yeo7_fsaverage5()

%--------------------------------------------------------------------------
% CBIG_MFMem_rfMRI_DrawOverlap_WI_Desikan68_Yeo7_fsaverage5
%
% show estimated W and I in Desikan68 parcellation of fsaverage5 with Yeo's 7-network boundary
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%--------------------------------------------------------------------------


%% read data (model parameter)
main_dir = pwd;
data_dir = fullfile(main_dir,'data'); 

cd('..')
high_dir = pwd;
lib_dir = fullfile(high_dir, 'lib');
addpath(lib_dir);
parameter_file_name = 'Example_Estimated_Parameter.mat';
Yeo7_label_file_name = '1000subjects_clusters007_ref.mat';
load ([data_dir '/' parameter_file_name],'Para_E');

Wvector = Para_E(1:68);
Ivector = Para_E(69:68*2);

colorname = 'cool';




%% load desikan68

lh_mesh_fsavg = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'aparc.annot');
rh_mesh_fsavg = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'aparc.annot');

lh_label_desikan = lh_mesh_fsavg.MARS_label;
rh_label_desikan = rh_mesh_fsavg.MARS_label;


%% load 7-network 

load ([data_dir '/' Yeo7_label_file_name]);

lh_label_desikan_data = lh_label_desikan;
rh_label_desikan_data = rh_label_desikan + 36;

 
lh_label_desikan_data(lh_label_desikan_data == 1) = 0;
lh_label_desikan_data(lh_label_desikan_data == 5) = 0;

rh_label_desikan_data(rh_label_desikan_data == 37) = 0;
rh_label_desikan_data(rh_label_desikan_data == 41) = 0;


%% Draw w, self connection
%create new color table
%generate a color table
HW = CBIG_MFMem_rfMRI_Num2Color(Wvector,colorname);
HW = squeeze(HW(:,1,:));
HW = HW.*255;
index = 1:1:72;
index([1,5,37,41]) = [];
HW_1 = ones(72,3);
HW_1(index,:) = HW;
HW_1 = [1 1 1; HW_1; 1 1 1];

%draw 
CBIG_DrawSurfaceMapsWithBoundary(lh_label_desikan_data,rh_label_desikan_data, lh_labels, rh_labels, 'fsaverage5', 'inflated',0,73, HW_1)
colorbar('off')
title('overlay W with 7-network boundary')

%generate a color bar
H_max = max(Wvector);
H_min = min(Wvector);
HW_colorbar = CBIG_MFMem_rfMRI_Num2Color(linspace(H_min,H_max,size(Wvector,1)),colorname);
figure
image(HW_colorbar)
set(gca,'xticklabel',{},'yticklabel',{})

%% draw I
HW = CBIG_MFMem_rfMRI_Num2Color(Ivector,colorname);
HW = squeeze(HW(:,1,:));
HW = HW.*255;
index = 1:1:72;
index([1,5,37,41]) = [];
HW_1 = ones(72,3);
HW_1(index,:) = HW;
HW_1 = [1 1 1; HW_1; 1 1 1];

%draw 
CBIG_DrawSurfaceMapsWithBoundary(lh_label_desikan_data,rh_label_desikan_data, lh_labels, rh_labels, 'fsaverage5', 'inflated',0,73, HW_1)
colorbar('off')
title('overlay I with 7-network boundary')

%generate a color bar
H_max = max(Ivector);
H_min = min(Ivector);
HW_colorbar = CBIG_MFMem_rfMRI_Num2Color(linspace(H_min,H_max,size(Wvector,1)),colorname);

figure
image(HW_colorbar)
set(gca,'xticklabel',{},'yticklabel',{})
rmpath(lib_dir);
end
