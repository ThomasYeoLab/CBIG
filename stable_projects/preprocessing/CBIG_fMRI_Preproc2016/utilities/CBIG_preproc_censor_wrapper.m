function CBIG_preproc_censor_wrapper(BOLD_in, outlier_file, TR, BOLD_interm_out, BOLD_final_out, loose_mask, max_mem, low_f, high_f)

% CBIG_preproc_censor_wrapper(BOLD_in, outlier_file, TR, BOLD_interm_out, BOLD_final_out, loose_mask, max_mem, low_f, high_f)
%
% Motion scrubbing for fMRI preprocessing. Users can perform bandpass
% filtering simultaneously by specifying low_f and high_f, where [low_f,
% high_f] (inclusive) is the passband. If low_f and high_f are not passed
% in, then bandpass filtering will not be performed.
% 
% This function uses the same censoring interpolation method as Power et
% al. 2014. Given the input fMRI file (BOLD_in), the outlier file
% (outlier_file), and TR, this function calls CBIG_preproc_censor.m to
% perform signal recovery and interpolation, where Lomb-Scargle Periodogram
% method is used. The intermediate signals and the final signals are stored
% in nifti files BOLD_interm_out and BOLD_final_out respectively.
%
% Note: At the beginning, this function detrend on the uncensored frames.
% This trend is not added back at the end. Moreover, no matter what option
% combination does the user pass in, the mean of the signal is not added
% back, which means if low_f = 0, 0 is NOT included in the passband.
%
% Input:
%     - BOLD_in: 
%       the BOLD file name before motion scrubbing (full path), e.g.
%       'subject_dir/subject_name/bold/subject_name_bld002_rest_skip4_stc_mc.nii.gz'
%
%     - outlier_file: 
%       file name of outliers (full path). In this file, each line is a
%       number of 0 or 1, where 0 indicates high motion frames that the
%       users want to censor, while 1 indicates low motion frames. e.g.
%       'subject_dir/subject_name/qc/subject_name_bld002_FDRMS0.2_DVARS50_motion_outliers.txt'
%
%     - TR: 
%       a string, the repetion time of fMRI data. The unit is milisecond.
%       e.g. '3000'
%
%     - BOLD_interm_out: 
%       the BOLD file name of intermediate result (full path). The
%       intermediate result means the signal of each time point (including
%       low motion time points) recovered by Lomb-Scargle Periodogram.
%       e.g. '<subject_dir>/<subject_name>/bold/<subject_name>_bld002_rest_skip4_stc_mc_interp_inter_FDRMS0.2_DVARS50.nii.gz'
%
%     - BOLD_out: 
%       the final BOLD file name after motion scrubbing (full path). 
%       (a) If bandpass filtering is not performed (low_f and high_f are not
%       passed in), the final signals are the interpolationed signal with
%       uncensored (low-motion) frames replaced by original signals. 
%       (b) If bandpass filtering is performed (both low_f and high_f are
%       passed in), the coefficients of frequency components outside the
%       passband are set to be 0, and then we use the masked coefficients
%       to recover the final signals.
%       e.g. '<subject_dir>/<subject_name>/bold/<subject_name>_bld002_rest_skip4_stc_mc_interp_FDRMS0.2_DVARS50.nii.gz'
% 
%     - loose_mask:
%       the filename of a loose whole brain mask (full path). The 
%       interpolation will only be done for the voxels within this loose
%       mask to save time (if passed in). The voxels outside the mask will
%       be set to 0.
%       e.g. '<subject_dir>/<subject_name>/bold/mask/<subject_name>.loosebrainmask.bin.nii.gz'
%       In the cases that you are not able to pass in a loose mask (for
%       instance, "BOLD_in" is a CIFTI dtseries file), you have two options:
%       (1) if you still want to pass in the latter two parameters "low_f"
%       and "high_f", you need to pass 'NONE' to "loose_mask" argument.
%       (2) if you do not need "low_f" and "high_f", you can skip
%       "loose_mask" argument as well.
% 
%     - max_mem:
%       a string of numbers to specify the maximal memory usage, or 'NONE'
%       (does not specify maximal memory usage). The unit is in Gigabyte.
%       In our code, we use a parameter k to adjust the maximal memory
%       usage, which is calculated according to this equation (concluded
%       from our tests)
%                max_mem (G) = 1 + (8e-4) * k * T
%       k is defined to determine the number of voxels processed each
%       time under this equation
%                V0 = floor(k * V / T / (oversample_fac/2))
%       where V0 is the number of voxels processed each time, V is the
%       total number of voxels within the whole brain (or within grey
%       matter, if "loose_mask" is passed in), T is the number of frames,
%       oversample_fac is an oversampling factor used in the Lomb-Scargle
%       algorithm, which is set to be 8.
%       We define this complicated equation is because we found the maximal
%       memory ussage is linearly proportional to the number of voxels
%       processed each time, and quadractically proportional to the number
%       of frames. 
%       We suggest the users to pass in a number that is 1G less than the
%       memory you will require from your job scheduler.
%       If 'NONE' is passed in, we use the default k = 20.
%
%     - low_f:
%       a string, low cut-off frequency. The passband includes cut-off
%       frequencies.
%       e.g. if the passband is [0, 0.08], then low_f is '0'.
%
%     - high_f:
%       a string, high cut-off frequency. The passband includes cut-off
%       frequencies. If the user wants to do highpass filtering, then
%       high_f is 'Inf'.
%       e.g. if the passband is [0, 0.08], then high_f is '0.08'.
%
% Example:
% CBIG_preproc_censor_wrapper('subject_dir/subject_name/bold/subject_name_bld002_rest_skip4_stc_mc.nii.gz',
%                             'subject_dir/subject_name/qc/subject_name_bld002_FDRMS0.2_DVARS50_motion_outliers.txt',
%                             '3000',
%                             'subject_dir/subject_name/bold/subject_name_bld002_rest_skip4_stc_mc_interp_inter_FDRMS0.2_DVARS50.nii.gz',
%                             'subject_dir/subject_name/bold/subject_name_bld002_rest_skip4_stc_mc_interp_FDRMS0.2_DVARS50.nii.gz',
%                             'subject_dir/subject_name/bold/mask/subject_name.loosebrainmask.nii.gz'
%                             '5'
%                             '0', '0.08')
% 
% Reference:
% 1) Power, Jonathan D., et al. "Methods to detect, characterize, and
%    remove motion artifact in resting state fMRI." Neuroimage 84 (2014):
%    320-341.
%
% Date: Jun.3, 2016
%
% Written by Jingwei Li.
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check loose mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_flag = 1;
if(~exist('loose_mask', 'var') || strcmp(loose_mask, 'NONE'))
    mask_flag = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check low_f and high_f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('low_f', 'var') && ~exist('high_f', 'var'))
    bandpass_flag = 0;
    fprintf('Do not perform bandpass filtering in censoring.\n');
elseif(exist('low_f', 'var') && exist('high_f', 'var'))
    bandpass_flag = 1;
    fprintf('Perform bandpass filtering in censoring.\n');
    fprintf('Passband is [%s, %s] (inclusive).\n', low_f, high_f);
else
    error('low_f or high_f does not exist! Please check the input arguments.\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read BOLD data and outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[input, in_vol, in_size] = read_fmri(BOLD_in);

% Read loose mask and apply it, if there is one.
if(mask_flag == 1)
    [~, mask_vol, ~] = read_fmri(loose_mask);
    mask_ind = find(mask_vol ~= 0);
    in_vol = in_vol(mask_ind, :);
end

outliers = dlmread(outlier_file);                    % N (number of timepoints) x 1 binary vector, 0 means censored (high-motion) frames.
outliers = ~outliers;                                % N x 1 binary vector, 0 means uncensored (low-motion) frames.

% detrend, trend is computed from uncensored frames
[in_vol, ~, ~, retrend] = CBIG_glm_regress_matrix(in_vol', [], 1, ~outliers);
in_vol = in_vol';

if(length(outliers)~=in_size(end))
    error('length of outlier file is not the same as the number of timepoints of input volume');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(in_vol, 2);

% check maximal memory usage
if(~exist('max_mem') || strcmp(max_mem, 'NONE'))
    k = 20;
else
    max_mem = str2num(max_mem);
    k = (max_mem - 1) / (8e-4) / N;                  % Tests show that 1 + (8e-4) * k * num_frames = max_mem (G)
end
fprintf('The factor used to split voxel batches is k = %f.\n', k);

TR = str2num(TR);
TR = TR/1000;
t = (1:N)' * TR;                                     % N x 1 vector, time of all frames
t_uncen = t(outliers==0);                            % (N-M) x 1 vector, time of uncensored frames

oversample_fac = 8;                                  % oversampling factor for Lomb-Scargle periodogram

if(bandpass_flag == 1)
    low_f = str2num(low_f);
    high_f = str2num(high_f);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lomb-Scargle Periodogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% divide voxels into branches, to reduce memory usage
voxbinsize = floor(k * size(in_vol, 1) / N / (oversample_fac/2));
voxbin = 1:voxbinsize:size(in_vol,1);
voxbin = [voxbin size(in_vol,1)+1];         % voxbin is the starting voxels in each branch

for v = 1:length(voxbin)-1
    fprintf('Dealing with voxels from %d to %d ...\n', voxbin(v), voxbin(v+1)-1);
    
    uncen_series = in_vol(voxbin(v):(voxbin(v+1)-1), outliers==0)';            % grab uncensored frames and the voxels in the voxel branch
    
    % interpolation, where interm_out_series is always the signal only with interppolation
    if(bandpass_flag == 0)
        % if no bandpass, out_series is the same signal with interm_out_series
        [out_vol(:, voxbin(v):(voxbin(v+1)-1)), interm_out_vol(:, voxbin(v):(voxbin(v+1)-1))] = CBIG_preproc_censor(uncen_series, t_uncen, t, oversample_fac, outliers);
    else
        % if bandpass, out_series is the interpolated signal within passband
        [out_vol(:, voxbin(v):(voxbin(v+1)-1)), interm_out_vol(:, voxbin(v):(voxbin(v+1)-1))] = CBIG_preproc_censor(uncen_series, t_uncen, t, oversample_fac, outliers, low_f, high_f);
    end
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interm_out_vol = single(interm_out_vol');
out_vol = single(out_vol');
if(bandpass_flag == 0)
    % if no bandpass, replace the uncensored frames with the original signal.
    out_vol(:, outliers==0) = in_vol(:, outliers==0);
end


% If there is a loose mask, construct output volumes
if(mask_flag == 1)
    tmp_interm_out = interm_out_vol;
    interm_out_vol = zeros(prod(in_size(1:3)), in_size(4));
    interm_out_vol(mask_ind, :) = tmp_interm_out;
    clear tmp_interm_out
    
    tmp_out = out_vol;
    out_vol = zeros(prod(in_size(1:3)), in_size(4));
    out_vol(mask_ind, :) = tmp_out;
    clear tmp_out
end

% write out output volumes
write_fmri(BOLD_interm_out, input, interm_out_vol, in_size);

write_fmri(BOLD_final_out, input, out_vol, in_size);



end





function [fmri, vol, vol_size] = read_fmri(fmri_name)

% [fmri, vol] = read_fmri(fmri_name)
% Given the name of functional MRI file (fmri_name), this function read in
% the fmri structure and the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%
% Output:
%     - fmri:
%       The structure read in by MRIread() or ft_read_cifti(). To save
%       the memory, fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) is
%       set to be empty after it is transfered to "vol".
%
%     - vol:
%       A num_voxels x num_timepoints matrix which is the content of
%       fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) after reshape.
%
%     - vol_size:
%       The size of fmri.vol (NIFTI) or fmri.dtseries (CIFTI).

if (isempty(strfind(fmri_name, '.dtseries.nii')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = single(fmri.vol);
    vol_size = size(vol);
    if(length(vol_size) < 4)
        vol = reshape(vol, prod(vol_size(1:3)), 1);
    else
        vol = reshape(vol, prod(vol_size(1:3)), vol_size(4));
    end
    fmri.vol = [];
else
    % if input file is CIFTI file
    fmri = ft_read_cifti(fmri_name);
    vol = single(fmri.dtseries);
    vol_size = size(vol);
    fmri.dtseries = [];
end


end



function write_fmri(fmri_name, fmri, vol, vol_size)

% function write_fmri(fmri_name, fmri, vol)
% This function write out a fmri strucure (fmri) with signal content (vol)
% into fmri_name.
% 
% Input:
%     - fmri_name:
%       The output fMRI file name (full path).
%
%     - fmri:
%       The structure for MRIwrite() or ft_write_cifti() to write out.
%
%     - vol:
%       The content of fMRI signals that need to be assgined to fmri.vol
%       (for NIFTI) or fmri.dtseries (for CIFTI).
%
%     - vol_size:
%       The size of volume when it was initially read in.

if(isempty(strfind(fmri_name, '.dtseries.nii')))
    % if output file is NIFTI file
    vol = reshape(vol, vol_size);
    fmri.vol = single(vol);
    MRIwrite(fmri, fmri_name);
else
    % if output file is CIFTI file
    vol = reshape(vol, vol_size);
    fmri.dtseries = single(vol);
    fmri_name = regexprep(fmri_name, '.dtseries.nii', '');
    ft_write_cifti(fmri_name, fmri, 'parameter', 'dtseries');
end


end
