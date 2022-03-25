function CBIG_preproc_filter_out_respiratory_pseudomotion(bold_file, ...
    motion_file, motion_regressors, out_FD, f_min, f_max)

% MR_filt = CBIG_preproc_motion_filtering(motion, TR, f_min, f_max)
%
% This function filters out respiratory pseudomotion from motion parameters
%
% Inputs:
%
%   - bold_file:
%     Path of the fMRI bold_file. We need this file to compute TR (Repetition time).
%
%   - motion_file:
%     A text file contains the N by 6 matrix of the motion parameters. N is
%     the number of frames. The first 3 columns should be rotation in radians, 
%     the last 3 columns should be translation in mm.
%
%   - motion_regressors
%     A text file contains the N by 6 matrix of the motion parameters. Both
%     this file and motion_file contains motion parameters by running
%     mcflirt on each run separately. The difference is that the motion
%     parameters in motion_file use the middle frame of each run as
%     reference frame. The motion parameters in this file use the first
%     frame of first run as reference.
%
%   - out_FD:
%     Path of the output FDRMS file after motion filtering. A text file
%     contrain the FDRMS values for each frame.
%
%   - f_min & f_max:
%     Both f_min and f_max should be string that can be converted to numbers or 
%     empty string. E.g., '0.3' or ''
%     f_min is the start freuqency of respiration anf f_max is the stop
%     frequency of the respiration. If f_max is empty, we use low pass filter 
%     and f_min is the stop frequency; if both f_min and f_max are provided
%     we use bandstop filter and the stopband is [f_min,f_max]
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%
% Reference: Fair et al., Correction of respiratory artifacts in MRI head motion parameters

%% input sanity check
TR = CBIG_preproc_infer_TR(bold_file);

if strcmp(f_min, '')
    error('f_min cannnot be empty');
end

%% load motion file and filter out respiratory pseudomotion
motion = load(motion_file);
motion_filt = CBIG_preproc_motion_filtering(motion, TR, f_min, f_max);

%% compute FDRMS from motion parameters and save out results
FDRMS = CBIG_preproc_compute_FDRMS_from_motion_parameters(motion_filt);
dlmwrite(out_FD, FDRMS, '\n');

%% filter out respiratory pseudomotion from motion regressors
motion = load(motion_regressors);
motion_filt = CBIG_preproc_motion_filtering(motion, TR, f_min, f_max);
copyfile(motion_regressors, [motion_regressors '.unfiltered']);
dlmwrite(motion_regressors, motion_filt, "  ");

end
