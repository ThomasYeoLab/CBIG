function colorscale = CBIG_MMP_GenerateColorscale(cmap, out_dir)

% function colorscale = CBIG_MMP_GenerateColorscale(cmap, out_dir)
% CBIG_MMP_create_annot_r_overlay_Schaefer400(input_vec, out_dir, lims, cmap)
% This function was was adapted from CBIG_TRBPC_create_annot_r_overlay_Schaefer400,
% creates the required colormap and saves the colorbar.
%
% Input:
% - out_dir: path to an output directory
% - cmap: a string that specifies what colour map is needed. Options are
% 'hot', 'cold', 'red_yell', 'blue_cyan' and 'hot_cold'.
%
% Output:
% colorscale: new colormap
%
% Written by Ruby Kong & CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

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
    if strcmp(cmap, 'blue_brown')
        RGB_new = [40 150 150;
            70  200 200;
            110 220 220;
            150 240 240;
            255 255 255;
            230 200 160;
            200 170 120;
            170 140 80;
            150 125 60]; 
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
    
    % print colorbar (WARNING: colorbar labels are hardcoded)
    % settings for horizontal
    hf = figure('Units','normalized'); 
    colormap(colorscale);
    cb = colorbar('north', 'XTickLabel',{'-1','-0.5','0','0.5','1'}, ...
               'XTick', [0,0.25,0.5,0.75,1]);
    set(gca,'Visible',false, 'fontsize', 20)
    cb.Position = [0.15 0.3 0.74 0.4];
    hf.Position(4) = 0.1000;
    saveas(gcf, fullfile(out_dir, 'colorbar_horz.png'))
    close(gcf)
    % settings for vertical
    hf = figure('Units','normalized'); 
    colormap(colorscale);
    cb = colorbar('west', 'XTickLabel',{'-1','-0.5','0','0.5','1'}, ...
               'XTick', [0,0.25,0.5,0.75,1],'YAxisLocation','right');
    set(gca,'Visible',false, 'fontsize', 20)
    cb.Position = [0.2 0.12 0.3 0.8];
    hf.Position(1) = 0;
    hf.Position(3) = 0.08;
    saveas(gcf, fullfile(out_dir, 'colorbar_vert.png'))
    close(gcf)
end