function command = CBIG_ViewRHAnnotationInFreeview(rh_annot_file, azimuth, elevation, output_file)
  % command = CBIG_ViewLHAnnotationInFreeview(rh_annot_file, azimuth, elevation, output_file)
  %
  % Produce a bash command to view a surface annotation of the right hemisphere in FreeView and capture the screenshot
  %
  % Input:
  %   - rh_annot_file  : .annot file of the right hemisphere to be visualized in FreeView
  %   - azimuth        : position of the FreeView's viewpoint about the z-axis. For the right hemisphere, setting azimuth = 180 shows the lateral view while azimuth = 0 shows the medial view
  %   - elevaton       : position of the FreeView's viewpoint about the y-axis. Setting elevation = 0 shows the side view of the brain, elevation = 90 shows the dorsal view, elevation = -90 shows the ventral view
  %   - output_file    : screenshot from FreeView
  % Output:
  %   - command        : bash command to open FreeView with the input annotation file with the view set by azimuth and elevation, and the screenshot saved as output_file
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  underlay = fullfile(getenv('FREESURFER_HOME'), 'subjects/fsaverage/surf/rh.inflated');
  command = ['freeview -f ' underlay ':annotation=' rh_annot_file ...
      ':edgethickness=0:color=lightgray' ...
      ' -viewport 3d' ...
      ' -cam azimuth ' num2str(azimuth) ' elevation ' num2str(elevation) ...
      ' -ss ' output_file ];
