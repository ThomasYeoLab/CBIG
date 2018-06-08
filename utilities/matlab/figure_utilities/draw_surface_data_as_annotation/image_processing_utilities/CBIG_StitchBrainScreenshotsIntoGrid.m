function command = CBIG_StitchBrainScreenshotsIntoGrid(lh_lateral_file, lh_medial_file, ...
    rh_lateral_file, rh_medial_file, dorsal_file, ventral_file, output_file)

% command = CBIG_StitchBrainScreenshotsIntoGrid(lh_lateral_file, lh_medial_file, ...
%     rh_lateral_file, rh_medial_file, dorsal_file, ventral_file, output_file)
%
% Combine screenshots of the brain into a grid layout
%
% Input:
%   - lh_lateral_file, lh_medial_file: screenshots of the lateral and medial views of the left hemisphere
%   - rh_lateral_file, rh_medial_file: screenshots of the lateral and medial views of the right hemisphere
%   - dorsal_file, ventral_flie      : screenshot of the dorsal view of both hemispheres
% Output:
%   - command                                    : Imagick command to perform the operation
%
% Example:
%   command = CBIG_StitchBrainScreenshotsIntoGrid('lh_lateral.png', 'lh_medial.png', ...
%      'rh_lateral.png', 'rh_medial.png', 'dorsal.png', 'ventral.png', 'grid.png')
%   Combine the input images into a grid layout as follows:
%   ---------------------------------------
%   | lh_lateral    rh_lateral   dorsal   |
%   | lh_medial     rh_medial    ventral  |
%   ---------------------------------------
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    command = ['convert -background white \( -background white ' lh_lateral_file ' ' rh_lateral_file ' ' dorsal_file ' +append \)' ...
        ' \( -background white ' lh_medial_file ' ' rh_medial_file ' ' ventral_file ' +append \) ' ...
        ' -append ' output_file];
