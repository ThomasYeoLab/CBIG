function command = CBIG_CombineImagesVertically(top_image, bottom_image, output_image)
% command = CBIG_CombineImagesVertically(top_image, bottom_image, output_image)
%
% Input:
%   - top_image, bottom_image: input images
%   - output_image           : final image with the two input images concatenated together vertically
% Output:
%   - command                : Imagick command to perform the concatenation
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  command = ['convert ' top_image ' ' bottom_image ' -background white -append ' output_image];
