function mni_coords = CBIG_RF_fsaverageVertex2MNICoord(hemi, vertices)
% mni_coords = CBIG_RF_fsaverageVertex2MNICoord(hemi, vertices)
%
% This function takes in a vertex number or a list of vertex numbers, to convert to the 
% corresponding RAS coordiantes in MNI152 space.
%
% Note that this conversion uses the MNI152-to-fsaverage mapping. Therefore the results 
% may not be consistent with using the final projection scripts to project fsaverage 
% data to MNI152 space.
%
% Input:
%     - hemi    :
%                hemisphere of the vertices ('lh' or 'rh')
%     - vertices:
%                the vertex number to convert, or an array of vertex
%                numbers to convert
%
% Output:
%     - coords:
%              3xN matrix containing the RAS coordinates in MNI152 space corresponding to 
%              the N input vertex numbers
%
% Example:
% mni_coords = CBIG_RF_fsaverageVertex2MNICoord('lh', 1:163842)
% This command finds the corresponding coordinates in MNI152 space for all the 1 to 163842 vertices 
% in the left hemisphere of fsaverage surface
% Note that you should be in the same directory as this script to run this command.
%
% Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Function usage
if nargin < 1
    disp('usage: mni_coords = CBIG_RF_fsaverageVertex2MNICoord(hemi, vertices)');
    return
end

% Load mappings
dir_uti = fileparts(fileparts(mfilename('fullpath')));
if strcmp(hemi, 'lh')
    load(fullfile(dir_uti, 'final_warps_FS5.3', 'lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat'), 'ras');
elseif strcmp(hemi, 'rh')
    load(fullfile(dir_uti, 'final_warps_FS5.3', 'rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat'), 'ras');
end

% Get corresponding coordinates
mni_coords = ras(:, vertices);

end
