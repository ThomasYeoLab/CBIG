function CBIG_ResizeImage(image_file, width, height, resized_file)

% CBIG_ResizeImage(image_file, width, height, resized_file)
%
% Resize an image to a new dimensions
%
% Input:
%   - image_file  : original image_file
%   - width       : width of the output image. If width = 0, scale the
%   image to the desired height and retains the original image ratio.
%   - height      : height of hte ouput image. If height = 0, scale the
%   image to the desired width and retains the original image ratio.
%   - resized_file: resized image file 
%
% Example:
%   CBIG_ResizeImage('lh_lateral.png', 1200, 500, 'lh_lateral.resized.png')
%   Resize the iamge 'lh_lateral.png' to the size 1200x500 pixel and save
%   the new image as 'lh_lateral.resized.png'
%   
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if width == 0
        system(['convert ' image_file ' -resize x' num2str(height) ' ' resized_file]);
    end
    
    if height == 0
        system(['convert ' image_file ' -resize ' num2str(width) ' ' resized_file]);
    end
    
    if width > 0 && height > 0
        system(['convert ' image_file ' -resize ' num2str(width) 'x' num2str(height) ' ' resized_file]);
    end
