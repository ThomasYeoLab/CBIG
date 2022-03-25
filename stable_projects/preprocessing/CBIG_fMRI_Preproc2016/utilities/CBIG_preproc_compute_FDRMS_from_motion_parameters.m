function FDRMS = CBIG_preproc_compute_FDRMS_from_motion_parameters(motion)

% FDRMS = CBIG_ABCD_proc_compute_FDRMSinson(motion)
% 
% This function compute the FDRMS (Jenkinson) from the 6 motion parameters
% 
% Inputs:
%   - motion:
%     #frame*6 matrix. The 6 motion parameters. 
%     The first 3 should be rotation in radians, the last 3 should be translation in mm.
%
% Oupouts:
%   - FDRMS:
%     #frame * 1 vector. reletive FDRMS for each frame
%
% Reference: Yan, Chao-Gan, et al. 
% "A comprehensive assessment of regional variation in the impact of head micromovements on functional connectomics." 
% Neuroimage 76 (2013): 183-201.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

T = cell(size(motion,1),1);
for i = 1:size(motion,1)
    alpha = motion(i,1);
    beta = motion(i,2);
    gamma = motion(i,3);
    R1 = [1 0 0 0
        0 cos(alpha) sin(alpha) 0
        0 -sin(alpha) cos(alpha) 0
        0 0 0 1];
    R2 = [cos(beta) 0 sin(beta) 0
        0 1 0 0
        -sin(beta) 0 cos(beta) 0
        0 0 0 1];
    R3 = [cos(gamma) sin(gamma) 0 0
        -sin(gamma) cos(gamma) 0 0
        0 0 1 0
        0 0 0 1];
    t = [1 0 0 motion(i,4)
        0 1 0 motion(i,5)
        0 0 1 motion(i,6)
        0 0 0 1];
    T{i} = t*R1*R2*R3;
end

FDRMS = zeros(length(T),1);
for i = 2:length(T)
    T_t = T{i}/(T{i-1}) - eye(4);
    A = T_t(1:3,1:3);
    b = T_t(1:3,4);
    FDRMS(i) = sqrt(1/5*80*80*trace(A*A')+b'*b);
end
        
end
