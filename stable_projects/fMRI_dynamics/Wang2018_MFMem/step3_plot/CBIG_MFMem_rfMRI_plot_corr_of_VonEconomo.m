function CBIG_MFMem_rfMRI_plot_corr_of_VonEconomo(save_dir_input)

%--------------------------------------------------------------------------
% CBIG_MFMem_rfMRI_plot_corr_of_VonEconomo()
%
% show the correlation of estimated W and I in Desikan68 parcellation with
% the Von Economo data. 
% W and I are averaged across the right and left hemishpere to martch the 
% Von Economo data.
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%--------------------------------------------------------------------------

%% input parameters
% save_dir_input: Output file directory

%% read data (model parameter)
main_dir = pwd;
data_dir = fullfile(main_dir,'data'); 

cd('..')
high_dir = pwd;
lib_dir = fullfile(high_dir, 'lib');
save_dir = save_dir_input;
addpath(lib_dir);
parameter_file_name = 'Example_Estimated_Parameter.mat';
SC_file_name = 'FCSC_Desikan68_Raphael_Wang.mat';
load ([data_dir '/' parameter_file_name],'Para_E');

%average cross right and left hemisphere
Wvector = Para_E(1:68);
Wvector34 = (Wvector(1:34)+Wvector(35:68))/2;
Ivector = Para_E(69:68*2);
Ivector34 = (Ivector(1:34)+Ivector(35:68))/2;
load([data_dir '/' SC_file_name],'SC_test');
SC = SC_test;
SCvector = sum(SC,2); %caculate all the incoming fibers (NOC: number of count)
SCvector34 = (SCvector(1:34)+SCvector(35:68))/2;


%% read von Economo data
load ([data_dir '/ECONOMO_data_Martijn.mat'])

Layer3_cellcontent = [layerIIItotal_cell_content_mm3_aparc'];
Layer3_cellsize = [layerIIItotal_cell_content_cellsize_aparc'];
Layer3_thickness = [layerIIItotal_thickness_overall_aparc];
Layer5_cellsize = [layerVtotal_cell_content_cellsize_aparc'];
Layer5_cellcontent = [layerVtotal_cell_content_mm3_aparc'];
Layer5_thickness = [layerVtotal_thickness_overall_aparc];
Layer1_cellsize = [layerItotal_cell_content_cellsize_aparc'];
Layer1_cellcontent = [layerItotal_cell_content_mm3_aparc'];
Layer1_thickness = [layerItotal_thickness_overall_aparc];
Layer6_cellsize = [layerVItotal_cell_content_cellsize_aparc'];
Layer6_cellcontent = [layerVItotal_cell_content_mm3_aparc'];
Layer6_thickness = [layerVItotal_thickness_overall_aparc];


Layer4_cellsize = [layerIVtotal_cell_content_cellsize_aparc'];
Layer4_cellcontent = [layerIVtotal_cell_content_mm3_aparc'];
Layer4_thickness = [layerIVtotal_thickness_overall_aparc];
Layer4_nan_indx = isnan(Layer4_cellsize);
Layer4_cellsize(isnan(Layer4_cellsize)) = 0;  %"only region2,5,15,23,25"
Layer4_cellcontent(isnan(Layer4_cellcontent)) = 0;%"only region2,5,15,23,25"
Layer4_thickness(isnan(Layer4_thickness)) = 0;

Layer2_cellsize = [layerIItotal_cell_content_cellsize_aparc'];
Layer2_cellcontent = [layerIItotal_cell_content_mm3_aparc'];
Layer2_thickness = [layerIItotal_thickness_overall_aparc];
Layer2_nan_indx = isnan(Layer2_cellsize);
Layer2_cellsize(isnan(Layer2_cellsize)) = 0;  %"only region2"
Layer2_cellcontent(isnan(Layer2_cellcontent)) = 0; %"only regin2"
Layer2_thickness(isnan(Layer2_thickness)) = 0;

layernum = 6*ones(34,1);
layernum = layernum-Layer2_nan_indx-Layer4_nan_indx;
LayerT_cellsize = (Layer1_cellsize+Layer2_cellsize+Layer3_cellsize+Layer4_cellsize+Layer5_cellsize+Layer6_cellsize)./layernum;
LayerT_cellcontent = (Layer1_cellcontent+Layer2_cellcontent+Layer3_cellcontent+Layer4_cellcontent+Layer5_cellcontent+Layer6_cellcontent)./layernum;


%% computing correlation

disp('pearson')
disp('        W  I  SC')
disp(['layer1   1  8 15'; 'layer2   2  9 16'; 'layer3   3 10 17'; 'layer4   4 11 18'; 'layer5   5 12 19'; 'layer6   6 13 20'; 'layerT   7 14 21'])


corr_matrix_full = cell(15,7);
corr_matrix_full{1,1} = 'Desikan68 R+L';
corr_matrix_full{2,1} = 'Layer1_cellDensity';
corr_matrix_full{3,1} = 'Layer2_cellDensity';
corr_matrix_full{4,1} = 'Layer3_cellDensity';
corr_matrix_full{5,1} = 'Layer4_cellDensity';
corr_matrix_full{6,1} = 'Layer5_cellDensity';
corr_matrix_full{7,1} = 'Layer6_cellDensity';
corr_matrix_full{8,1} = 'LayerT_cellDensity';
corr_matrix_full{9,1} = 'Layer1_cellSize';
corr_matrix_full{10,1} = 'Layer2_cellSize';
corr_matrix_full{11,1} = 'Layer3_cellSize';
corr_matrix_full{12,1} = 'Layer4_cellSize';
corr_matrix_full{13,1} = 'Layer5_cellSize';
corr_matrix_full{14,1} = 'Layer6_cellSize';
corr_matrix_full{15,1} = 'LayerT_cellSize';
corr_matrix_full{1,2} = 'W';
corr_matrix_full{1,3} = 'P-value';
corr_matrix_full{1,4} = 'I';
corr_matrix_full{1,5} = 'P-value';
corr_matrix_full{1,6} = 'SC_NOS';
corr_matrix_full{1,7} = 'P-value';

[r_value,p_value] = corrcoef(Wvector34(Layer1_cellcontent>0),Layer1_cellcontent(Layer1_cellcontent>0));
corr_matrix_full{2,2} = r_value(1,2);
corr_matrix_full{2,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer2_cellcontent>0),Layer2_cellcontent(Layer2_cellcontent>0));
corr_matrix_full{3,2} = r_value(1,2);
corr_matrix_full{3,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer3_cellcontent>0),Layer3_cellcontent(Layer3_cellcontent>0));
corr_matrix_full{4,2} = r_value(1,2);
corr_matrix_full{4,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer4_cellcontent>0),Layer4_cellcontent(Layer4_cellcontent>0));
corr_matrix_full{5,2} = r_value(1,2);
corr_matrix_full{5,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer5_cellcontent>0),Layer5_cellcontent(Layer5_cellcontent>0));
corr_matrix_full{6,2} = r_value(1,2);
corr_matrix_full{6,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer6_cellcontent>0),Layer6_cellcontent(Layer6_cellcontent>0));
corr_matrix_full{7,2} = r_value(1,2);
corr_matrix_full{7,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(LayerT_cellcontent>0),LayerT_cellcontent(LayerT_cellcontent>0));
corr_matrix_full{8,2} = r_value(1,2);
corr_matrix_full{8,3} = p_value(1,2);

[r_value,p_value] = corrcoef(Wvector34(Layer1_cellsize>0),Layer1_cellsize(Layer1_cellsize>0));
corr_matrix_full{9,2} = r_value(1,2);
corr_matrix_full{9,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer2_cellsize>0),Layer2_cellsize(Layer2_cellsize>0));
corr_matrix_full{10,2} = r_value(1,2);
corr_matrix_full{10,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer3_cellsize>0),Layer3_cellsize(Layer3_cellsize>0));
corr_matrix_full{11,2} = r_value(1,2);
corr_matrix_full{11,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer4_cellsize>0),Layer4_cellsize(Layer4_cellsize>0));
corr_matrix_full{12,2} = r_value(1,2);
corr_matrix_full{12,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer5_cellsize>0),Layer5_cellsize(Layer5_cellsize>0));
corr_matrix_full{13,2} = r_value(1,2);
corr_matrix_full{13,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(Layer6_cellsize>0),Layer6_cellsize(Layer6_cellsize>0));
corr_matrix_full{14,2} = r_value(1,2);
corr_matrix_full{14,3} = p_value(1,2);
[r_value,p_value] = corrcoef(Wvector34(LayerT_cellsize>0),LayerT_cellsize(LayerT_cellsize>0));
corr_matrix_full{15,2} = r_value(1,2);
corr_matrix_full{15,3} = p_value(1,2);

[r_value,p_value] = corrcoef(Ivector34(Layer1_cellcontent>0),Layer1_cellcontent(Layer1_cellcontent>0));
corr_matrix_full{2,4} = r_value(1,2);
corr_matrix_full{2,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer2_cellcontent>0),Layer2_cellcontent(Layer2_cellcontent>0));
corr_matrix_full{3,4} = r_value(1,2);
corr_matrix_full{3,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer3_cellcontent>0),Layer3_cellcontent(Layer3_cellcontent>0));
corr_matrix_full{4,4} = r_value(1,2);
corr_matrix_full{4,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer4_cellcontent>0),Layer4_cellcontent(Layer4_cellcontent>0));
corr_matrix_full{5,4} = r_value(1,2);
corr_matrix_full{5,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer5_cellcontent>0),Layer5_cellcontent(Layer5_cellcontent>0));
corr_matrix_full{6,4} = r_value(1,2);
corr_matrix_full{6,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer6_cellcontent>0),Layer6_cellcontent(Layer6_cellcontent>0));
corr_matrix_full{7,4} = r_value(1,2);
corr_matrix_full{7,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(LayerT_cellcontent>0),LayerT_cellcontent(LayerT_cellcontent>0));
corr_matrix_full{8,4} = r_value(1,2);
corr_matrix_full{8,5} = p_value(1,2);

[r_value,p_value] = corrcoef(Ivector34(Layer1_cellsize>0),Layer1_cellsize(Layer1_cellsize>0));
corr_matrix_full{9,4} = r_value(1,2);
corr_matrix_full{9,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer2_cellsize>0),Layer2_cellsize(Layer2_cellsize>0));
corr_matrix_full{10,4} = r_value(1,2);
corr_matrix_full{10,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer3_cellsize>0),Layer3_cellsize(Layer3_cellsize>0));
corr_matrix_full{11,4} = r_value(1,2);
corr_matrix_full{11,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer4_cellsize>0),Layer4_cellsize(Layer4_cellsize>0));
corr_matrix_full{12,4} = r_value(1,2);
corr_matrix_full{12,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer5_cellsize>0),Layer5_cellsize(Layer5_cellsize>0));
corr_matrix_full{13,4} = r_value(1,2);
corr_matrix_full{13,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(Layer6_cellsize>0),Layer6_cellsize(Layer6_cellsize>0));
corr_matrix_full{14,4} = r_value(1,2);
corr_matrix_full{14,5} = p_value(1,2);
[r_value,p_value] = corrcoef(Ivector34(LayerT_cellsize>0),LayerT_cellsize(LayerT_cellsize>0));
corr_matrix_full{15,4} = r_value(1,2);
corr_matrix_full{15,5} = p_value(1,2);

[r_value,p_value] = corrcoef(SCvector34(Layer1_cellcontent>0),Layer1_cellcontent(Layer1_cellcontent>0));
corr_matrix_full{2,6} = r_value(1,2);
corr_matrix_full{2,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer2_cellcontent>0),Layer2_cellcontent(Layer2_cellcontent>0));
corr_matrix_full{3,6} = r_value(1,2);
corr_matrix_full{3,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer3_cellcontent>0),Layer3_cellcontent(Layer3_cellcontent>0));
corr_matrix_full{4,6} = r_value(1,2);
corr_matrix_full{4,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer4_cellcontent>0),Layer4_cellcontent(Layer4_cellcontent>0));
corr_matrix_full{5,6} = r_value(1,2);
corr_matrix_full{5,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer5_cellcontent>0),Layer5_cellcontent(Layer5_cellcontent>0));
corr_matrix_full{6,6} = r_value(1,2);
corr_matrix_full{6,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer6_cellcontent>0),Layer6_cellcontent(Layer6_cellcontent>0));
corr_matrix_full{7,6} = r_value(1,2);
corr_matrix_full{7,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(LayerT_cellcontent>0),LayerT_cellcontent(LayerT_cellcontent>0));
corr_matrix_full{8,6} = r_value(1,2);
corr_matrix_full{8,7} = p_value(1,2);

[r_value,p_value] = corrcoef(SCvector34(Layer1_cellsize>0),Layer1_cellsize(Layer1_cellsize>0));
corr_matrix_full{9,6} = r_value(1,2);
corr_matrix_full{9,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer2_cellsize>0),Layer2_cellsize(Layer2_cellsize>0));
corr_matrix_full{10,6} = r_value(1,2);
corr_matrix_full{10,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer3_cellsize>0),Layer3_cellsize(Layer3_cellsize>0));
corr_matrix_full{11,6} = r_value(1,2);
corr_matrix_full{11,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer4_cellsize>0),Layer4_cellsize(Layer4_cellsize>0));
corr_matrix_full{12,6} = r_value(1,2);
corr_matrix_full{12,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer5_cellsize>0),Layer5_cellsize(Layer5_cellsize>0));
corr_matrix_full{13,6} = r_value(1,2);
corr_matrix_full{13,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(Layer6_cellsize>0),Layer6_cellsize(Layer6_cellsize>0));
corr_matrix_full{14,6} = r_value(1,2);
corr_matrix_full{14,7} = p_value(1,2);
[r_value,p_value] = corrcoef(SCvector34(LayerT_cellsize>0),LayerT_cellsize(LayerT_cellsize>0));
corr_matrix_full{15,6} = r_value(1,2);
corr_matrix_full{15,7} = p_value(1,2);




%% plot

%draw celldensity
writeSize = 12;

figure
x = Wvector34(LayerT_cellcontent>0);
y = LayerT_cellcontent(LayerT_cellcontent>0);
scatter(x,y,'filled')
title('W and layerT neuronal density')
xlabel('Self-recurrent W')
ylabel('neuronal density')
p = polyfit(x,y,1);
yfit = p(1)*x+p(2);
hold on
plot(x,yfit,'r--');
set(gca,'fontsize',writeSize,'fontname','arial','TickDir','out');
saveas(gca,[save_dir '/W_cellDensity'],'epsc')
saveas(gca,[save_dir '/W_cellDensity'],'png')

figure
x = Ivector34(LayerT_cellcontent>0);
y = LayerT_cellcontent(LayerT_cellcontent>0);
scatter(x,y,'filled')
title('I and layerT neuronal density')
xlabel('Subcortical input I')
ylabel('Neuronal density')
p = polyfit(x,y,1);
yfit = p(1)*x+p(2);
hold on
plot(x,yfit,'r--');
set(gca,'fontsize',writeSize,'fontname','arial','TickDir','out');
saveas(gca,[save_dir '/I_cellDensity'],'epsc')
saveas(gca,[save_dir '/I_cellDensity'],'png')


%draw cellsize 

figure
x = Wvector34(LayerT_cellsize>0);
y = LayerT_cellsize(LayerT_cellsize>0);
scatter(x,y,'filled')
title('W and layerT neuronal size')
xlabel('Self-recurrent W')
ylabel('Neuronal size')
p = polyfit(x,y,1);
yfit = p(1)*x+p(2);
hold on
plot(x,yfit,'r--');
set(gca,'fontsize',writeSize,'fontname','arial','TickDir','out');
saveas(gca,[save_dir '/W_cellSize'],'epsc')
saveas(gca,[save_dir '/W_cellSize'],'png')


figure
x = Ivector34(LayerT_cellsize>0);
y = LayerT_cellsize(LayerT_cellsize>0);
scatter(x,y,'filled')
title('I and layerT neuronal size')
xlabel('Subcortical input I')
ylabel('Neuroanl size')
p = polyfit(x,y,1);
yfit = p(1)*x+p(2);
hold on
plot(x,yfit,'r--');
set(gca,'fontsize',writeSize,'fontname','arial','TickDir','out');
saveas(gca,[save_dir '/I_cellSize'],'epsc')
saveas(gca,[save_dir '/I_cellSize'],'png')



%% FDR test
%p-value inclusive SC
p_value_all = corr_matrix_full(2:15,[3 5 7]);
p_value_all = cell2mat(p_value_all);
p1 = p_value_all(:);
[significant_index_SC, final_threshold_SC] = FDR(p1, 0.05);
disp('inclusive Total layer:')
disp(['final threshold:' num2str(final_threshold_SC)]);
disp(['discovery count:' num2str(length(significant_index_SC))])
disp('discovery:')
significant_index_SC

%p-value exclusive SC
p2 = p_value_all(:,1:2);
[significant_index_noSC, final_threshold_noSC] = FDR(p2, 0.05);
disp('inclusive Total layer:')
disp('exclusive SC')
disp(['final threshold:' num2str(final_threshold_noSC)])
disp(['discovery count:' num2str(length(significant_index_noSC))])
disp('discovery:')
significant_index_noSC



%% save to xlsx file
wfilename = 'VonEconomo_correlation.xlsx';
xlswrite([save_dir '/' wfilename],corr_matrix_full,1,'A1');

wfilename = 'VonEconomo_correlation_desikan';
save([save_dir '/' wfilename],'corr_matrix_full','final_threshold_noSC','final_threshold_SC');
rmpath(lib_dir);
end










