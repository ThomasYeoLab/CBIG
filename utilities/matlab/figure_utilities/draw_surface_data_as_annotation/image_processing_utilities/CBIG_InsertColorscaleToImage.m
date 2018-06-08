function command = CBIG_InsertColorscaleToImage(colorscale_file, underlay_image, output_image, scale_ratio, x_offset_percentage, y_offset_percentage)

% command = CBIG_InsertColorscaleToImage(colorscale_file, underlay_image, output_image, scale_ratio, x_offset_percentage, y_offset_percentage)
%
% Produce a bash command to insert a colorscale to an image
%
% Input:
%   - colorscale_file: path to the image file of the colorscale
%   - underlay_iamge : path to the underlaying image for the colorscale to be added to
%   - output_image   : path to the resulting image
%   - scale_ratio    : the ratio to which the unerlaying image is scaled to. Default: 0.5
%   - x_offset_ratio : horizontal position of the colorscale as a percentage of the underlaying image's width. Default: 0.31
%   - y_offset_ratio : vertical position of the colorscale as a percentage of the underlaying image's height. Default: 0.49
% Output:
%   - command        : Imagick's script to insert the colorscale to the underlaying image
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if nargin < 4  
        scale_ratio = 0.5;
        x_offset_percentage = 0.31;
        y_offset_percentage = 0.49;
    end
    
    if isnumeric(scale_ratio)
        scale_ratio = num2str(scale_ratio);
    end
    
    if isnumeric(x_offset_percentage)
        x_offset_percentage = num2str(x_offset_percentage);
    end
    
    if isnumeric(y_offset_percentage)
        y_offset_percentage = num2str(y_offset_percentage);
    end
    
    command = ['w=`identify ' colorscale_file '| cut -f 3 -d " " | sed s/x.*//`; ' ...
        'new_w=`awk "BEGIN {printf(\"%.0f\", $w * ' scale_ratio ')}"`; ' ...
        'underlay_w=`identify ' underlay_image '| cut -f 3 -d " " | sed s/x.*//`; ' ...
        'colorscale_x=`awk "BEGIN {printf(\"%.0f\", $underlay_w * ' x_offset_percentage ')}"`; ' ...
        'h=`identify ' colorscale_file '| cut -f 3 -d " " | sed s/.*x//`; ' ...
        'new_h=`awk "BEGIN {printf(\"%.0f\", $h * ' scale_ratio ')}"`; ' ...
        'underlay_h=`identify ' underlay_image '| cut -f 3 -d " " | sed s/.*x//`; ' ...
        'colorscale_y=`awk "BEGIN {printf(\"%.0f\", $underlay_h * ' y_offset_percentage ')}"`; ' ...
        'composite -compose atop -geometry +${colorscale_x}+${colorscale_y} ' colorscale_file ' ' underlay_image ' ' output_image];  
