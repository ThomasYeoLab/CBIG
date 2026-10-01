function FD = CBIG_preproc_compute_FD_from_motion_parameters(motion)

% FD = CBIG_preproc_compute_FD_from_motion_parameters(motion)
% 
% This function computes FD (Power et al, 2012) from the 6 motion parameters.
% 
% Inputs:
%   - motion:
%     #frame*6 matrix. The 6 motion parameters. 
%     The first 3 columns should be rotation in radians,
%     the last 3 columns should be translation in mm.
%
% Outputs:
%   - FD:
%     #frame * 1 vector. relative FD for each frame
%
% Example:
%     motion = rand(N, 6); % random motion matrix (#frame x 6)
%     FD = CBIG_preproc_compute_FD_from_motion_parameters(motion)
%
% Reference: Power, et al.
% "Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion."
% Neuroimage 59 (2012): 2142-2154.
%
% Written by Trevor Tan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% compute translation contribution to FD
if size(motion, 2) ~= 6
    error('Input argument ''motion'' should have 6 columns (#frame x 6).');
end
trans = motion(:, 4:end);
trans_diff = abs(trans(1:end-1, :) - trans(2:end, :));
trans_sum = sum(trans_diff, 2);

% compute rotation contribution to FD
rot = motion(:, 1:3);
rot_diff = abs(rot(1:end-1, :) - rot(2:end, :));
rot_sum = sum(rot_diff, 2);
% multiply rotation (radians) by 50 mm as in Power et al, 2012
rot_sum = rot_sum .* 50;

% compute FD
FD = trans_sum + rot_sum;
% default set 1st frame FD = 0 mm
FD = [0; FD];

% number of frames sanity check
assert(length(FD) == length(motion), ...
    'ERROR: FD length (%d) not equal to number of frames (%d).', ...
    length(FD), length(motion));
end
