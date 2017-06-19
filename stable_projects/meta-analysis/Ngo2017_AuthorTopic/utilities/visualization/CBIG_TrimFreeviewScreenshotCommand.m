function command = TrimFreeviewScreenshotCommand(input_file, output_file, width_shaving_pct, height_shaving_pct)
  if nargin < 3
    width_shaving_pct = 0.30;
    height_shaving_pct = 0.20;
  end
  
  if isnumeric(width_shaving_pct)
    width_shaving_pct = num2str(width_shaving_pct);
  end
  
  if isnumeric(height_shaving_pct)
    height_shaving_pct = num2str(height_shaving_pct);
  end

  command = ['w=`identify ' input_file '| cut -f 3 -d " " | sed s/x.*//`; ' ...
      'shaving_w=`awk "BEGIN {printf(\"%.0f\", $w * ' width_shaving_pct ')}"`; ' ...
      'h=`identify ' input_file '| cut -f 3 -d " " | sed s/.*x//`; ' ...
      'shaving_h=`awk "BEGIN {printf(\"%.0f\", $h * ' height_shaving_pct ')}"`; ' ...
      'convert ' input_file ' -shave ${shaving_w}x${shaving_h} ' output_file];