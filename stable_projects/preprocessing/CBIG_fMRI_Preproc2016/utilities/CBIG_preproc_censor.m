function [out_series, interm_out_series, f, cos_coeff, sin_coeff] = CBIG_preproc_censor(in_series, in_sample_time,...
    out_sample_time, oversample_fac, outliers, low_f, high_f)

% [out_series, interm_out_series, f, cos_coeff, sin_coeff] =
% CBIG_preproc_censor(in_series, in_sample_time, out_sample_time, oversample_fac, low_f, high_f)
%
% This function estimates lost data point from unevenly sampled input data.
% Users can perform bandpass filtering simultaneously by specifying low_f
% and high_f, where [low_f, high_f] (inclusive) is the passband. If low_f
% and high_f are not passed in, then bandpass filtering will not be
% performed.
%
% This function computes estimated signal interm_out_series (at time
% out_sample_time) from unevenly sampled signal in_series (sampled at time
% in_sample_time) using Lomb-Scarge periodogram spectral analysis technique.
% The length of in_sample_time must be equal to the number of rows of
% in_series. The routine will calculate the coefficients for an increasing
% sequence of frequencies up to the Nyquist frequency, with an oversampling
% factor of loversample_fac (typically >= 4). 
% 
% If low_f and high_f are passed in, the coefficients of the frequencies
% outside the passband are set to be 0. out_series is recovered by the
% masked coefficients.
% 
% If low_f and high_f are not passed in, out_series is the same as
% interm_out_series.
%
% Input:
%     - in_series: 
%       Input signal matrix (treat each column as a channel) with dimension
%       of num_in_time x nChannel, where num_in_time is the number of
%       unevenly sampled time points, nChannel could be number of voxels in
%       fMRI data.
%
%     - in_sample_time: 
%       Nonuniform sampling time sequence with dimension of num_in_time x
%       1. in_series(k) is sampled at time in_sample_time(k), where k is
%       from 1 to num_in_time.
%
%     - out_sample_time: 
%       A vector, recovering time sequence with dimension of num_out_time x
%       1. The recovered signal out_series(k) is sampled at
%       out_sample_time(k), where k is from 1 to num_out_time.
%
%     - oversample_fac: 
%       A scalar, the oversampling factor (typically >= 4). The interval of
%       frequencies is 1/(T*oversample_fac), where T is the whole time
%       interval (T = max(in_sample_time) - min(in_sample_time)). e.g. 8
%
%     - outliers:
%       A num_out_time x 1 binary vector. 0 means low-motion frames; 1
%       means high-motion frames. To be clearer, the indices of 0s
%       correspond to in_sample_time.
%
%     - low_f:
%       A string, low cut-off frequency. The passband includes cut-off
%       frequencies. If the user does not want to perform bandpass
%       filtering, please do not pass in this parameter.
%       e.g. if the passband is [0, 0.08], then low_f is '0'.
%
%     - high_f:
%       A string, high cut-off frequency. The passband includes cut-off
%       frequencies. If the user wants to do highpass filtering, then
%       high_f is 'Inf'. If the user does not want to perform bandpass
%       filtering, please do not pass in this parameter.
%       e.g. if the passband is [0, 0.08], then high_f is '0.08'.
%
% Output:
%     - out_series:
%       Depend on whether bandpass filtering is performed. If bandpass
%       filtering is performed, the coefficients of frequencies outside the
%       passband are set to be 0. And out_series is recovered by the masked
%       coefficients. If bandpass filtering is not performed, out_series is
%       the same as interm_out_series. 
%       It has the dimension of num_out_time x nChannel.
%
%     - interm_out_series: 
%       Recovered signal matrix with dimension of num_out_time x nChannel.
%       The interpolation is done with the coefficients on all frequency
%       components.
%
%     - f: 
%       A num_freq_bin x 1 vector. It is an increasing sequence of
%       frequencies, and num_freq_bin is the number of frequencies in this
%       sequence. It is the frequency sequence that is considered when this
%       function tries to recover the signal.
%
%     - cos_coeff: 
%       A num_freq_bin x 1 x nChannel matrix. It is the coeffients of
%       cos(2*pi*f*in_sample_time - tau) components, where tau is the phase
%       delay defined in Lomb-Scargle periodogram. It is used to recover
%       signal.
%
%     - sin_coeff: 
%       A num_freq_bin x 1 x nChannel matrix. It is the coeffients of
%       sin(2*pi*f*in_sample_time - tau) components, where tau is the phase
%       delay defined in Lomb-Scargle periodogram. It is used to recover
%       signal.
% 
% Example:
% [out_sample_series, f, cos_coeff, sin_coeff] =
% CBIG_preproc_censor(in_sample_series, in_sample_time, out_sample_time, 8, outliers, 0, 0.08)
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

%% Check input variables
if size(in_sample_time,2) > 1
    error('Input argument ''in_sample_time'' should be a column vector');
end

if size(out_sample_time,2) > 1
    error('Input argument ''out_sample_time'' should be a column vector');
end

if size(outliers,2) > 1
    error('Input argument ''outliers'' should be a column vector');
end

%% check low_f and high_f
if(~exist('low_f', 'var') && ~exist('high_f', 'var'))
    bandpass_flag = 0;
elseif(exist('low_f', 'var') && exist('high_f', 'var'))
    bandpass_flag = 1;
else
    error('low_f or high_f does not exist! Please check the input arguments.\n');
end

%% compute dimension parameters
in_series = single(in_series);
num_in_time = length(in_sample_time);
nChannel = size(in_series,2);

%% compute mean and std of input series
input_mean = single(mean(in_series,1));                   % 1 x nChannel
input_std = single(std(in_series,1));                     % 1 x nChannel
in_series = bsxfun(@minus, in_series, input_mean);


%% frequency sequence
% T is the time period of all samples in uncensored frames.
% sample_period = T / num_in_time
% sample_frequency = 1 / sample_period
% sample_frequency / 2 is the uppper bound of frequency sequency in recovery
T = max(in_sample_time) - min(in_sample_time);
f = (1/(T*oversample_fac):1/(T*oversample_fac):num_in_time/(2*T)).';       % num_freq_bin x 1
w = single(2*pi*f);

%% sin and cos components
% num_freq_bin x 1
tau = atan2(sum(sin(2*w*in_sample_time.'),2), sum(cos(2*w*in_sample_time.'),2)) ./ (2*w); 
% num_freq_bin x num_in_time
cterm = single(cos(bsxfun(@times, w, bsxfun(@minus, in_sample_time.', tau))));  
sterm = single(sin(bsxfun(@times, w, bsxfun(@minus, in_sample_time.', tau))));

in_series = reshape(in_series, 1, num_in_time, nChannel);          % 1 x num_in_time x nChannel

%% coefficients of cos components (using equation (3) in supplementary material in the same folder)
cos_h = bsxfun(@times, in_series, cterm);                 % num_freq_bin x num_in_time x nChannel
numerator = sum(cos_h, 2);                                % num_freq_bin x 1 x nChannel
denominator = sum(cterm.^2, 2);                           % num_freq_bin x 1
cos_coeff = bsxfun(@times, numerator, 1./denominator);    % num_freq_bin x 1 x nChannel
clear numerator denominator cos_h cterm

%% coefficients of sin components (using equation (3) in supplementary material in the same folder)
sin_h = bsxfun(@times, in_series, sterm);
numerator = sum(sin_h, 2); 
denominator = sum(sterm.^2, 2);
sin_coeff = bsxfun(@times, numerator, 1./denominator);
clear numerator denominator sin_h sterm

%% recovery
cterm_new = single(cos(bsxfun(@times, w, bsxfun(@minus, out_sample_time.', tau))));
sterm_new = single(sin(bsxfun(@times, w, bsxfun(@minus, out_sample_time.', tau))));

% first, recover the intermediate result (no bandpass)
% use equation (4) in supplementary material in the same folder
% num_freq_bin x num_out_time x nChannel
interm_out_series = single(bsxfun(@times, cos_coeff, cterm_new) + bsxfun(@times, sin_coeff, sterm_new));  
% 1 x num_out_time x nChannel
interm_out_series = sum(interm_out_series,1);                     
% num_out_time x nChannel
interm_out_series = reshape(interm_out_series, length(out_sample_time), nChannel);  

output_std = single(std(interm_out_series(outliers==0, :), 1));            % 1 x nChannel
stdfac = input_std ./ output_std;
interm_out_series = bsxfun(@times, interm_out_series, stdfac);

% second, compute the final output (with bandpass or with replacement, depends on bandpass_flag)
if(bandpass_flag == 1)
    % if bandpass, recover the signal within passband
    % mask coefficients
    fmask = ((f < low_f) | (f > high_f));
    cos_coeff(fmask==1, 1, :) = 0;
    sin_coeff(fmask==1, 1, :) = 0;
    
    % use equation (4) in supplementary material in the same folder
    % num_freq_bin x num_out_time x nChannel
    out_series = single(bsxfun(@times, cos_coeff, cterm_new) + bsxfun(@times, sin_coeff, sterm_new));     
    % 1 x num_out_time x nChannel
    out_series = sum(out_series,1);                
    % num_out_time x nChannel
    out_series = reshape(out_series, length(out_sample_time), nChannel); 
    
    output_std = single(std(out_series(outliers==0, :), 1));                       % 1 x nChannel
    stdfac = input_std ./ output_std;
    out_series = bsxfun(@times, out_series, stdfac);
else
    % if no bandpass, out_series is the same with interm_out_series
    out_series = interm_out_series;
end


end

