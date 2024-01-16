function CBIG_MMP_create_annot_r_overlay_Schaefer400(input_vec, out_dir, lims, cmap)

% CBIG_MMP_create_annot_r_overlay_Schaefer400(input_vec, out_dir, lims, cmap)
% This function was was adapted from CBIG_TRBPC_create_annot_r_overlay_Schaefer400,
% projects a 400 length vector onto the Schaefer 400 parcellation. 
% Requires 'Schaefer2018_400Parcels_17Networks_fslr_164K.mat' to be in the same folder.
%
% Input:
% - input_vec: a 1x400 vector
% - out_dir: path to an output directory
% - lims: min and max values for colour scale e.g. [0 100]
% - cmap: a string that specifies what colour map is needed. Options are
% 'hot', 'cold', 'red_yell', 'blue_cyan' and 'hot_cold'.
%
% Written by Ruby Kong & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%SET min and max threshold for colorbar:
min_th = lims(1);
max_th = lims(2);

%SET annot file output dir:
mkdir(out_dir);

%SET annot file output filename:
lh_annot_file=fullfile(out_dir,'lh_data.annot');
rh_annot_file=fullfile(out_dir,'rh_data.annot');

%load Alex parcellation
load('Schaefer2018_400Parcels_17Networks_fslr_164K.mat');

lh_data=zeros(163842,1);
rh_data=zeros(163842,1);
input_vec(isnan(input_vec)) = 0;
for i=1:400
    lh_data(HCP.lh_labels==i)=input_vec(i);
    rh_data(HCP.rh_labels==i)=input_vec(i);
end

colorscale = CBIG_MMP_GenerateColorscale(cmap, out_dir);

[lh_convert_labels, lh_colortable] = CBIG_MMP_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(lh_data,...
    140, colorscale, min_th, max_th);
[rh_convert_labels, rh_colortable] = CBIG_MMP_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(rh_data,...
    140, colorscale, min_th, max_th);


%inset boundary color as black
boundary_rgb = [10 10 10];
boundary_rgb_sum = 10 + 10*2^8 + 10*2^16;
boundary_color = [boundary_rgb 0 boundary_rgb_sum];
lh_colortable.table = [lh_colortable.table; boundary_color];
rh_colortable.table = [rh_colortable.table; boundary_color];
lh_colortable.numEntries = 142;
rh_colortable.numEntries = 142;
lh_colortable.struct_names(142) = {'boundary'};
rh_colortable.struct_names(142) = {'boundary'};

lh_convert_labels(lh_data==0) = lh_colortable.table(end-1, 5);
rh_convert_labels(rh_data==0) = rh_colortable.table(end-1, 5);

lh_convert_labels(HCP.lh_labels==0) = lh_colortable.table(end, 5);
rh_convert_labels(HCP.rh_labels==0) = rh_colortable.table(end, 5);


write_annotation(lh_annot_file, 0:(length(lh_convert_labels)-1), lh_convert_labels, lh_colortable);
write_annotation(rh_annot_file, 0:(length(rh_convert_labels)-1), rh_convert_labels, rh_colortable);

end
