function []=CBIG_gwMRF_create_FSL_LUT(filename,colors)

% [] = CBIG_gwMRF_create_FSL_LUT(filename,colors)
%
% This function takes a nx3 colors vector, where n is the number of colors
% and creates a FSL compatible LUT file.

%input 
% - filename: the output .lut file
% - colors: nx3 colors vector, where n is the number of colors
% output is generated as files
%
%example
% CBIG_gwMRF_create_FSL_LUT([base_folder,'/MNI/Schaefer2018_',num2str(j),'Parcels_',num2str(i) ,'Networks_order.lut'],[lh_s.table(2:end,[1:3]);rh_s.table(2:end,[1:3])])
% Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

fileID=fopen(filename,'w');
fprintf(fileID,'%!VEST-LUT\n');
fprintf(fileID,'%%BeginInstance\n');
fprintf(fileID,'<<\n');
fprintf(fileID,'/SavedInstanceClassName /ClassLUT\n');
fprintf(fileID,'/PseudoColorMinimum 0.00\n');
fprintf(fileID,'/PseudoColorMaximum 1.00\n');
fprintf(fileID,'/PseudoColorMinControl /Low\n');
fprintf(fileID,'/PseudoColorMaxControl /High\n');

fprintf(fileID,'/PseudoColormap [\n');
for i=1:max(size(colors))
    fprintf(fileID,'<-color{%1.6f,%1.6f,%1.6f}->\n',colors(i,1)/255,colors(i,2)/255,colors(i,3)/255);
end
fprintf(fileID,']\n');
fprintf(fileID,'>>\n');
fprintf(fileID,'%%EndInstance\n');
fprintf(fileID,'%%EOF\n');

fclose(fileID);
