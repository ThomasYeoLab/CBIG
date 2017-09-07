function command = CBIG_PadImage(image_file, padded_file, x_padding, y_padding)
  % command = CBIG_PadImage(image_file, x_padding, y_padding)
  %
  % Produce a command to add transparent paddings to an image
  %
  % Input
  %   - image_file : image file to be padded
  %   - x_padding  : padding on the left and right borders
  %   - y_padding  : padding on the top and bottom borders
  %   - padded_file: padded image
  % Output:
  %   - command    : Imagick command to perform the padding
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if nargin < 3
    x_padding = 5;
    y_padding = 5;
  end
  
  command = ['convert ' image_file ' -matte -bordercolor white -border ' num2str(x_padding) 'x' num2str(y_padding) ' ' padded_file];
