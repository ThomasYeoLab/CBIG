function [min_w, min_h] = CBIG_GetMinimumDimensionsOfImages(image_file_list)
% [min_w, min_h] = CBIG_GetMinimumDimensionsOfImages(image_file_list)
%
% Get the minimum width and height from a series of images
%
% Input:
%   - image_file_list: list of input images
% Output:
%   - min_w          : the smallest width among the input images
%   - min_h          : the smallest height among the input images
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  w_list = zeros(length(image_file_list), 1);
  h_list = zeros(length(image_file_list), 1);
  
  count = 1;
  for image_file = image_file_list
    f = cell2mat(image_file);

    info = imfinfo(f);
    w = info.Width;
    h = info.Height;

    w_list(count) = w;
    h_list(count) = h;
    count = count + 1;
  end
  
  min_w = min(w_list);
  min_h = min(h_list);
