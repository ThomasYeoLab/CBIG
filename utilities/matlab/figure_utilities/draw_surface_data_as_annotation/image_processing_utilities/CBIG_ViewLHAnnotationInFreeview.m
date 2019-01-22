function command = CBIG_ViewLHAnnotationInFreeview(lh_annot_file, underlay_file, azimuth, elevation, output_file)

% command = CBIG_ViewLHAnnotationInFreeview(lh_annot_file, underlay_file, azimuth, elevation, output_file)
%
% Produce a bash command to view a surface annotation of the left hemisphere in FreeView and capture the screenshot
%
% Input:
%   - lh_annot_file  : .annot file of the left hemisphere to be visualized in FreeView
%   - underlay_file  : file containing the underlaying brain surface
%   - azimuth        : position of the FreeView's viewpoint about the z-axis.
%                      For the left hemisphere, setting azimuth = 0 shows the lateral
%                      view while azimuth = 180 shows the medial view
%   - elevaton       : position of the FreeView's viewpoint about the y-axis.
%                      Setting elevation = 0 shows the side view of the brain,
%                      elevation = 90 shows the dorsal view, elevation = -90
%                      shows the ventral view
%   - output_file    : screenshot from FreeView
% Output:
%   - command        : bash command to open FreeView with the input annotation file
%                      with the view set by azimuth and elevation, and the screenshot
%                      saved as output_file
%
% Example:
%   lh_underlay_file = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage5', 'surf', 'lh.inflated');
%   command = CBIG_ViewLHAnnotationInFreeview(lh_component.annot, lh_underlay_file, 180, 0, 'lh_medial.png')
%   Return bash command to save the screenshot of the medial view of the
%   left hemisphere. The surface annotation is saved in
%   'lh_component.annot'. The underlay surface is the inflated fsaverage5
%   surface of the left hemisphere. The medial view is shown in FreeView by
%   setting the azimuth to 180 and elevation to 0. The output screenshot is
%   saved to 'lh_medial.png'
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    command = ['freeview -f ' underlay_file ':annotation=' lh_annot_file ...
        ':edgethickness=0' ...
        ' -viewport 3d' ...
        ' -cam azimuth ' num2str(azimuth) ' elevation ' num2str(elevation) ...
        ' -ss ' output_file ];
