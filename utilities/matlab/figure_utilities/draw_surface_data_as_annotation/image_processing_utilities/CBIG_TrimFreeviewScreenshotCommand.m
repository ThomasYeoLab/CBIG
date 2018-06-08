function command = CBIG_TrimFreeviewScreenshotCommand(input_file, output_file)

% command = CBIG_TrimFreeviewScreenshotCommand(input_file, output_file)
%
% Produce a bash command to trim transparent paddings surrounding an image
%
% Input:
%   - input_file : original image with the transparent padding to be removed
%   - output_file: image wiht the padding removed
% Output:
%   - command    : Imagick command to perform the trimming
%
% Example:
%   command = CBIG_TrimFreeviewScreenshotCommand('lh_lateral.png', 'lh_lateral.trimmed.png')
%   Return command to trim the transparent background at the borders of the image
%   lh_lateral.png and save the new image as lh_lateral_trimmed.png
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    command = ['convert ' input_file ' -trim +repage ' output_file];
