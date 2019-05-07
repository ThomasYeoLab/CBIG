function colorscale = CBIG_GenerateClearbrainColorscale(res, min_val, max_val, output_dir)

% colorscale = CBIG_GenerateClearbrainColorscale(res, min_val, max_val, output_dir)
%
% Produce a discretized colorscale based on Human Connectome Workbench's
% clear_brain color pallete
% (https://www.humanconnectome.org/software/workbench-command/-all-commands-help).
%
% Input:
%   - res       : number of discrete values in the colorscale
%   - min_val   : minimum of the colorscale
%   - max_val   : maximum of the colorscale
%   - output_dir: output directory of the colorscale's image
% Output:
%   - colorscale: colorscale in Matlab's format with res number of values
%   discretized from the clear_brain color scheme 
%
% colorscale = CBIG_GenerateClearbrainColorScale(28, 5e-5, 1, ...
%   '/Work/data/colorscale')
%   Generate a discrete colorscale of 28 bins with the value ranging from
%   5e-5 to 1. The colors come from the Human Connectome Workbench's
%   clear_brain color pallete.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    RGB_COLORSCALE = [ 255 0 0;
                       255 105 0;
                       255 153 0;
                       255 255 0;
                       16  176 16;
                       0 255 0;
                       127 127 204;
                       76 76 127;
                       51 51 76;
                       102 0 51 ];
    
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
    
    Y_OFFSET = 30;
    X_OFFSET = 100;
    WIDTH = 3;
    HEIGHT = 10;
    for i = 1:res
        c = colorscale(i, :);
        rectangle('Position', [X_OFFSET + i * WIDTH, Y_OFFSET, WIDTH, HEIGHT], 'FaceColor', c, 'EdgeColor', c);
    end
    axis([0, size(colorscale, 1) * WIDTH + X_OFFSET * 2, 0, Y_OFFSET * 2 + HEIGHT]);
    set(gca, 'Visible', 'off');
    
    min_label = num2str(min_val);
    max_label = num2str(max_val);
    label_y = round(Y_OFFSET + 0.4 * HEIGHT);
    label_font_size = 25;
    text(X_OFFSET - 60, label_y, min_label, 'FontSize', label_font_size);
    text(X_OFFSET + size(colorscale, 1) * WIDTH + 8, label_y, max_label, 'FontSize', label_font_size);
    
    if nargin == 4
        set(gcf, 'PaperPositionMode', 'auto');
        print(fullfile(output_dir, 'colorscale'), '-dpng', '-r0');
    end
