function motion_filt = CBIG_preproc_motion_filtering(motion, TR, f_min, f_max)

% MR_filt = CBIG_preproc_motion_filtering(motion, TR, f_min, f_max)
%
% This function filters out respiratory pseudomotion from motion estimates
%
% Inputs:
%
%   - motion:
%     N by 6 matrix of the motion estimates. N is the number of frames
%
%   - TR:
%     TR of the fMRI image in seconds.
%
%   - f_max & f_min:
%     If f_max is empty, we use low pass filter and f_min is the stop frequency;
%     if both f_min and f_max are provided we use bandstop filter and the
%     stopband is [f_min,f_max]
%
% Outputs:
%
%   - motion_filt
%     N by 6 filtered motion estimates. N is the number of frames
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%
% Reference: Fair et al., Correction of respiratory artifacts in MRI head motion estimates

%% input sanity check
order = 4;
f_min = str2num(f_min);
f_max = str2num(f_max);

if isempty(f_min)
    error('f_min (start frequency of respiration) cannot be empty');
end

if isempty(f_max)
    filt_type = 'lp';
else
    filt_type = 'notch';
end

if ~isempty(f_max) && f_min >= f_max
    error('f_min should be less than f_max')
end

%% filter design
fs = 1/TR;
fNy = fs/2;

switch filt_type

    case 'lp'

        Wn = f_min/fNy;
        b_filt = fir1(order, Wn, 'low');
        a_filt = 1;
        num_f_apply = 0;

    case 'notch'

        W_notch = [f_min,f_max]/fNy;
        Wn = mean(W_notch);
        bw = diff(W_notch);
        fprintf("%s, %s", Wn, bw)
        [b_filt, a_filt] = iirnotch(Wn, bw);
        num_f_apply = floor(order / 2); % if order<4 apply filter 1x, if order=4 2x, if order=6 3x

end

%% Read individual movement regressors files
motion_filt = filtfilt(b_filt,a_filt,motion);
for i=1:num_f_apply-1
    motion_filt = filtfilt(b_filt,a_filt,motion_filt);
end

end