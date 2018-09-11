function CBIG_MFMem_rfMRI_DrawBarchart_W_Desikan68_into_Yeo7_fsaverage5()

%-------------------------------------------------------------------------
% CBIG_MFMem_rfMRI_DrawBarchart_W_Desikan68_into_Yeo7_fsaverage5
% 
% show mean value of W in Yeo's 7-network
% what to do?
% 1) give the value of W to each voxel, according to the Desikan labels, mapping W from Desikan68 to voxel
% 2) caculate the value for Yeo 52 components, mapping W from each voxel to Yeo 52 components
% 3) caculate the value for Yeo 7-networks, mapping W from Yeo 52 components to 7-networks
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%--------------------------------------------------------------------------


%% read data set
%read data (model parameter)
main_dir = pwd;
data_dir = fullfile(main_dir,'data'); 
parameter_file_name = 'Example_Estimated_Parameter.mat';
Thomas52_label_file_name = 'ThomasYeo_52_components_read_out.mat';
Yeo7_label_file_name = '1000subjects_clusters007_ref.mat';
load ([data_dir '/' parameter_file_name],'Para_E');

Wvector = Para_E(1:68);
Ivector = Para_E(69:68*2);


%% read labels, Thomas Yeo 52 componets
load ([data_dir '/' Thomas52_label_file_name]);

lh_labels_fs = zeros(10242,1);
rh_labels_fs = zeros(10242,1);

lh_labels_fs(lh_labels_fs_thomas52(:,1)) = lh_labels_fs_thomas52(:,2);
rh_labels_fs(rh_labels_fs_thomas52(:,1)) = rh_labels_fs_thomas52(:,2);
rh_labels_fs(rh_labels_fs>0) = rh_labels_fs(rh_labels_fs>0) + 26;


%% read labels, find labels(ThomasYeo07; Deskian) for each voxel
%read 7-networks
load ([data_dir '/' Yeo7_label_file_name]);

%%  read desikan
lh_mesh_fsavg = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'aparc.annot');
rh_mesh_fsavg = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'aparc.annot');

lh_label_desikan = lh_mesh_fsavg.MARS_label';
rh_label_desikan = rh_mesh_fsavg.MARS_label';

lh_label_desikan_data = lh_label_desikan;
rh_label_desikan_data = rh_label_desikan + 36;

lh_label_desikan_data(lh_label_desikan_data == 1) = 0;
lh_label_desikan_data(lh_label_desikan_data == 5) = 0;

rh_label_desikan_data(rh_label_desikan_data == 37) = 0;
rh_label_desikan_data(rh_label_desikan_data == 41) = 0;


%% give the value of W to each voxel, according to the Desikan labels
%generate new lh_Ivector_36, rh_Ivector_36
lh_Wvector_36 = zeros(36,1);
rh_Wvector_36 = zeros(36,1);
lh_Wvector_36([2:4 6:36]) = Wvector(1:34);  %Wvector has no No.1 unknown, and No.5 corpus callsum
rh_Wvector_36([2:4 6:36]) = Wvector(35:68);

%give the value of W to each voxel, according to the Desikan
lh_vox_value = zeros(10242,1);
rh_vox_value = zeros(10242,1);
for i = 1:10242
    lh_vox_value(i) = lh_Wvector_36(lh_label_desikan(i));
end
for i = 1:10242
    rh_vox_value(i) = rh_Wvector_36(rh_label_desikan(i));
end


%% caculate the value for for Thomas Yeo 52 components

for i = 1:(length(unique([lh_labels_fs rh_labels_fs]))-1)
   
    lh_vox_indx_for_52components{i} = find(lh_labels_fs == i);
    rh_vox_indx_for_52components{i} = find(rh_labels_fs == i);
   
    
    lh_vox_value_tmp = lh_vox_value(lh_vox_indx_for_52components{i});
    lh_vox_value_tmp = lh_vox_value_tmp(lh_vox_value_tmp>0);
    lh_vox_num_for_52components = length(lh_vox_value_tmp);
    
    rh_vox_value_tmp = rh_vox_value(rh_vox_indx_for_52components{i});
    rh_vox_value_tmp = rh_vox_value_tmp(rh_vox_value_tmp>0);
    rh_vox_num_for_52components = length(rh_vox_value_tmp);
        
    lr_value_52{i} = [lh_vox_value_tmp;rh_vox_value_tmp];

    lr_value_52components(i) = (sum(lh_vox_value_tmp)+sum(rh_vox_value_tmp))/(lh_vox_num_for_52components+rh_vox_num_for_52components);
    
    lr_size_52components(i) =  lh_vox_num_for_52components+rh_vox_num_for_52components;
end


%% assign 52 components to 7 networks
%read components name
Name52 = [lh_labels_name; rh_labels_name];
% assign 17 network name
Name7 = {'Vis';'SomMot';'DorsAttn';'SalventAttn';'Limbic';'Cont';'Default'};
% assign 114 components to 17 network
label_52_to_7 = zeros(51,1);
label_52_to_7(1) = 1;
label_52_to_7(2) = 2;
label_52_to_7(3:5) = 3;
label_52_to_7(6:10) = 4;
label_52_to_7(11:12) = 5;
label_52_to_7(13:21) = 6;
label_52_to_7(22:26) = 7;
label_52_to_7(27) = 1;
label_52_to_7(28) = 2;
label_52_to_7(29:31) = 3;
label_52_to_7(32:37) = 4;
label_52_to_7(38:39) = 5;
label_52_to_7(40:46) = 6;
label_52_to_7(47:51) = 7;


%% caculate W mean and error bar in 7-networks 
for i = 1:7

        regions_network7{i} = find(label_52_to_7 == i);
        regions_network7_value{i} = lr_value_52components(regions_network7{i});
        regions_network7_size{i} = lr_size_52components(regions_network7{i});
        regions_network7_value_mean(i) = mean(regions_network7_value{i});
        regions_network7_value_std(i) = std(regions_network7_value{i});
        regions_network7_errorbar(i) = regions_network7_value_std(i)/sqrt(length(regions_network7{i}));
        
end


%% plot statistic

figure
bar(regions_network7_value_mean,'FaceColor','w');
ylim([0.4 0.8])
hold on
errorbar([1:7], regions_network7_value_mean, regions_network7_errorbar,'+k');
set(gca,'xTick',[1:7],'xTicklabel', Name7, 'yTick',[0.4:0.1:0.7])
set(gca,'Fontsize',12,'Fontname','arial','TickDir','out')
hold off
title('Self-recurrent W')

end
