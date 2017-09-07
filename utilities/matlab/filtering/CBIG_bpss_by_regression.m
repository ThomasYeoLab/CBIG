function output_matrix = CBIG_bpss_by_regression(input_matrix, low_f, high_f, sample_period, censor)

% output_matrix = CBIG_bpss_by_regression(input_matrix, low_f, high_f, sample_period, censor)
%
% Using linear regression to perform bandpass filtering.
%
% Given the low frequency (low_f) and high frequency (high_f) of the
% passband, this function first construct a frequency mask where 1 denotes
% the frequency to be removed, whereas 0 denotes the frequency that is
% kept. According to the frequency mask, this funtion construct cosine and
% sine regressors correponding to the frequencies to be removed and call
% function CBIG_glm_regress_matrix.m to regress out the components of these
% frequencies in input_matrix. The residual is returned to be the bandpass
% filtered signal, which only contains the components of frequencies that
% the users are interested into. The sample_period is used to construct the
% frequency sequence. The boundary frequencies low_f and high_f are not
% thrown.
%
% This function can also handle a censoring vector (censor) with the same
% length of timeseries in input_matrix. The element of this vector is
% either 0 or 1, where 1 means this time point has low motion and the
% frequencies of this time point will be considered, while 0 means this
% frame has high motion so that the user wants to exclude it.
%
% NOTE: low_f < high_f. Bandstop is not allowed in this function.
%
% Input:
%     - input_matrix:
%       a T x N data matrix, where T is number of time points, N is number
%       of samples. input_matrix(t, n) means the signal intensity of n-th
%       sample in time t.
%
%     - low_f:
%       a scalar, low cutoff frequency (Hz). e.g. 0.001.
%       
%     - high_f:
%       a scalar, high cutoff frequency (Hz). e.g. 0.08.
% 
%     - sample_period:
%       a scalar, the sample period (second). e.g. 2.
% 
% Output:
%     - output_matrix:
%       a T x N matrix, where T is number of time points, N is number of
%       samples. output_matrix(t, n) means the signal intensity of n-th
%       sample in time t. output_matrix is the bandpass filtered results.     
%
% Example:
% output_matrix = CBIG_bandpass_by_regression(input_data, 0.001, 0.08, 2, censor)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



% check if low_f < high_f
if(low_f >= high_f)
    error('ERROR: Bandstop is not permitted.')
end

% check if high pass frequency is positive
if(high_f <= 0)
    error('ERROR: high pass frequency must be positive.')
end

% check if "censor" in passed in
if(~exist('censor', 'var'))
    censor = [];
end

% sample frequency
Fs = 1/sample_period;

% sample points
L = size(input_matrix, 1);
nf = floor(L/2);

if (mod(L, 2) == 0)
    % If L is even
    f = [0:(L/2-1) L/2:-1:1]*Fs/L;  %frequency vector
elseif (mod(L, 2) == 1)
    % If L is odd
    f = [0:(L-1)/2 (L-1)/2:-1:1]*Fs/L;
end

% construct frequency mask
% 1 means the corresponding frequency needs to be filtered out
% 0 means keeping this frequency
fmask = ((f < low_f) | (f > high_f));
fmask = fmask(1:nf+1);
demean = -1;
if(fmask(1))
    demean = 0;
end
fmask(1) = 0;

% construct regressors
pp = 1:nf;
pp = pp(fmask(2:end)==1);                      % the indices of frequencies that to be removed
fq = 2 * pi * pp / L;                        
rads = bsxfun(@times, [0:L-1]', fq);
cosf = cos(rads);
sinf = sin(rads);
if(mod(L,2) == 0 && fmask(end) == 1)
    sinf = sinf(:, 1:end-1);
end

% regress out the stopband frequencies
input_matrix = double(input_matrix);
output_matrix = CBIG_glm_regress_matrix(input_matrix, [cosf sinf], demean, censor);

end
