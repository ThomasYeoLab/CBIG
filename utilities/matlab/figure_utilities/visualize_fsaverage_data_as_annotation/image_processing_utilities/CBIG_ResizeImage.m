function CBIG_ResizeImage(image_file, width, height, resized_file)
  % CBIG_ResizeImage(image_file, width, height)
  %
  % Resize an image to a new dimensions
  %
  % Input:
  %   - image_file  : original image_file
  %   - width       : width of the output image
  %   - height      : height of hte ouput image
  %   - resized_file: resized image file 
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
