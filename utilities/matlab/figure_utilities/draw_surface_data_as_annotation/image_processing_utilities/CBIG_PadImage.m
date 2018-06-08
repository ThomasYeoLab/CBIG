function command = CBIG_PadImage(image_file, padded_file, x_padding, y_padding)

% command = CBIG_PadImage(image_file, padded_file, x_padding, y_padding)
%
%
% Produce a command to add transparent paddings to an image
%
% Input
%   - image_file : image file to be padded
%   - x_padding  : padding in pixels on the left and right borders. Default = 5
%   - y_padding  : padding in pixels on the top and bottom borders. Default = 5
%   - padded_file: padded image
% Output:
%   - command    : Imagick command to perform the padding
%
% Example:
%   command = CBIG_PadImage('lh_lateral.png', 'padded_lh_lateral.png', ... 
%      10, 10)
%   Return command 'convert lh_lateral.png -matte -bordercolor white -border
%     10x10 padded_lh_lateral.png' that adds a padidng of 10 pixels around
%     the image lh_lateral.png and saved the new image in 'padded_lh_lateral.png'
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if nargin < 3
        x_padding = 5;
        y_padding = 5;
    end
    
    command = ['convert ' image_file ' -matte -bordercolor white -border ' num2str(x_padding) 'x' num2str(y_padding) ' ' padded_file];
