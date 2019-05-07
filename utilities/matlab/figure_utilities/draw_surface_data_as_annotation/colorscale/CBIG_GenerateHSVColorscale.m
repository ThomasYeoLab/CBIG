function colorscale = CBIG_GenerateHSVColorscale(res, min_val, max_val, output_dir)

% colorscale = CBIG_GenerateHSVColorscale(res, min_val, max_val, output_dir)
%
% Produce a discretized colorscale based on Matlab's HSV color scheme
%
% Input:
%   - res       : number of discrete values in the colorscale
%   - min_val   : minimum of the colorscale
%   - max_val   : maximum of the colorscale
%   - output_dir: output directory of the colorscale's image
% Output:
%   - colorscale: colorscale in Matlab's format with res number of values discretized from the HSV color scheme 
%
% colorscale = CBIG_GenerateHSVColorscale(28, 5e-5, 1, ...
%   '/Work/data/colorscale')
%   Generate a discrete colorscale of 28 bins with the value ranging from
%   5e-5 to 1. The colors come from the Matlab's HSV color pallete.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    padding = round(res / 4);
    colorscale = colormap(hsv(res + padding));
    colorscale = flipud(colorscale(1:res, :));
    shade = 1.0;
    colorscale = colorscale * shade;
    
    Y_OFFSET = 30;
    X_OFFSET = 200;
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
    text(40, label_y, min_label, 'FontSize', label_font_size);
    text(X_OFFSET + size(colorscale, 1) * WIDTH + 30, label_y, max_label, 'FontSize', label_font_size);
    
    if nargin == 4
        set(gcf, 'PaperPositionMode', 'auto');
        print(fullfile(output_dir, 'colorscale'), '-dpng', '-r0');
    end
    
    close all;
