function colortable = CBIG_GenerateAnnotationColortable(discretization_res, colorscale)
  % colortable = CBIG_GenerateAnnotationColortable(discretization_res, colorscale)
  %
  % Generate a colortable suitable for surface annotation based on an input colorscale
  %
  % Input:
  %   - discretization_res: number of distinct discrete colors in the output color table
  %   - colorscale        : the input color scheme that the output colortable for the annotation is produced
  % Output:
  %   - colortable        : the output colortable with the format suitable for surface annotation
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % gray color of the underlaying brain
  gray_rgb = [200 200 200];
  gray_rgb_sum = 200 + 200*2^8 + 200*2^16;
  underlay_color = [gray_rgb 0 gray_rgb_sum];

  colortable.numEntries = discretization_res+1;
  colortable.orig_tab = '';
  colortable.struct_names = cell(discretization_res+1, 1);
  for i = 1:discretization_res+1
    colortable.struct_names{i} = 'unknown';
  end
  
  table = zeros(discretization_res+1, 5);
  for i = 1:discretization_res
    table(i, 1:3) = uint8(colorscale(i, :) * 255 + 0.5);
    table(i, 5)   = table(i, 1) + table(i, 2) * 2^8 + table(i, 3) *2^16;
  end
  table(discretization_res+1, :) = underlay_color;
  colortable.table = table;
