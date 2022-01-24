function CBIG_TRBPC_create_annot_r_overlay_Schaefer400(input_vec, out_dir, lims, cmap)

%% This function will project a vector onto the Schaefer 400 parcellation

% Required inputs
% - input_vec: a 1x400 vector
% - out_dir: path to an output directory
% - lims: min and max values for colour scale e.g. [0 100]
% - cmap: a string that specifies what colour map is needed. Options are
% 'hot' or 'cold'
%
% Written by Ruby Kong & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%%

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

colorscale = CBIG_GenerateColorscale(cmap);

[lh_convert_labels, lh_colortable] = CBIG_TRBPC_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(lh_data,...
    140, colorscale, min_th, max_th);
[rh_convert_labels, rh_colortable] = CBIG_TRBPC_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(rh_data,...
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


function colorscale = CBIG_GenerateColorscale(cmap)
    res = 140;
    if strcmp(cmap, 'hot')
        RGB_new = [255 0 0;
            255 31  31;
            255 56  56;
            255 92  92;
            255 120 120;
            255 158 158;
            255 186 186;
            255 214 214;
            255 255 255];
    end
    if strcmp(cmap, 'cold')
        RGB_new = [0 0 255;
            31  31  255;
            56  56  255;
            92  92  255;
            120 120 255;
            158 158 255;
            186 186 255;
            214 214 255;
            255 255 255];
    end
    if strcmp(cmap, 'red_yell')
        RGB_new = [255 255 255;
            255 236 94;
            232 208 28;
            181 148 0;
            181 97  0;
            181 51  0;
            181 0   0;
            117 0   0;
            0   0   0];
    end
    if strcmp(cmap, 'blue_cyan')
        RGB_new = [255 255 255;
            179 255 251;
            43  224 215;
            43  182 224;
            43  155 224;
            43  118 224;
            13  75  184;
            0   43  117;
            0   0   0]; 
    end
    if strcmp(cmap, 'hot_cold')
        RGB_new = [250 241 182;
            250 208 70;
            255 140 0;
            255 0   0;
            0   0   0;
            0   0   255;
            0   140 255;
            70  208 250;
            182 241 250];
    end
            
    RGB_COLORSCALE = RGB_new;
    orig_colorscale_length = size(RGB_COLORSCALE, 1);
    colorscale = [];
    step = round(res / orig_colorscale_length);
    
    count = 1;
    for i = 1:orig_colorscale_length-1
        for j = 0:step-1
            colorscale(count, :) = RGB_COLORSCALE(i, :) * (step - j) + RGB_COLORSCALE(i+1, :) * j;
            colorscale(count, :) = colorscale(count, :) * 1.0 / step;
            count = count + 1;
        end
    end
    colorscale = [colorscale; RGB_COLORSCALE(orig_colorscale_length, :)];
    
    pad = (res - size(colorscale, 1)) / 2;
    for i = 1:pad
        colorscale = [RGB_COLORSCALE(1, :); colorscale];
    end
    
    pad = res - size(colorscale, 1);
    for i = 1:pad
        colorscale = [colorscale; RGB_COLORSCALE(orig_colorscale_length, :)];
    end
    
    colorscale = flipud(colorscale) / 255;
end