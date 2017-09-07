function command = CBIG_ReplaceBlackByWhiteBkgCommand(input_file, output_file)
% command = ReplaceBlackByWhiteBkgCommand(input_file, output_file)
%
% Replace the black background of FreeView screenshots with white color
%
% Input:
%   - inpute_file: original input image
%   - output_file: final image
% Output:
%   - command: Imagick command to clean up the image
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  command = ['convert ' input_file ...
      ' -alpha set -channel RGBA -fuzz 1% -fill white -floodfill +0+0 black ' ...
      output_file ];
