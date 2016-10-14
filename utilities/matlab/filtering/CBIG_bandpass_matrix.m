function output_mtx = CBIG_bandpass_matrix(input_mtx, low_f, high_f, sample_period)

% output_mtx = CBIG_bandpass_matrix(input_mtx, low_f, high_f, sample_period)
% 
% Using FFT to transform data (input_mtx) into frequency domain and set the 
% low frequency (low_f) and high frequency (high_f) to construct a bandpass 
% window [low_f, high_f], low_f and high_f frequency are included in this 
% window. Then this window will be applied to the data in the frequency 
% domain and use IFFT to transform the data back to temporal domain.
% In this function the sample period is set by sample_period.
%
% NOTE: low_f < high_f, because bandstop is not allowed in this function.
% 
% Input:
%     - input_mtx: 
%       a T x N data matrix, where T is number of time points, N is number of
%       samples. input_mtx(t, n) means the signal intensity of n-th sample 
%       in time t.
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
%     - output_mtx:
%       a T x N matrix, where T is number of time points, N is number of
%       samples. output_mtx(t, n) means the signal intensity of n-th sample 
%       in time t. output_mtx is the bandpass filtered results.     
%
% Example:
% output_mtx = CBIG_bandpass_matrix(input_data, 0.001, 0.08, 2)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% check whether low cutoff frequency is smaller that high frequency
% bandstop is not allowed
if (low_f >= high_f)
    error('ERROR: Bandstop is not allowed.')
end

% sample frequency
Fs = 1/sample_period;

% sample points
L = size(input_mtx, 1);


if (mod(L, 2) == 0)
    % If L is even
    f = [0:(L/2-1) L/2:-1:1]*Fs/L;  %frequency vector
elseif (mod(L, 2) == 1)
    % If L is odd
    f = [0:(L-1)/2 (L-1)/2:-1:1]*Fs/L;
end

% index of frequency that you want to keep
ind = ((low_f <= f) & (f <= high_f));

% create a rectangle window according to the cutoff frequency
rectangle = zeros(L, 1);
rectangle(ind) = 1;

% use fft to transform the time course into frequency domain
fprintf('FFT each time course.\n')
input_mtx_fft = fft(input_mtx);

% apply the rectangle window and ifft the signal from frequency domain to time domain
fprintf('Apply rectangle window and IFFT.\n')
output_mtx = ifft(bsxfun(@times, input_mtx_fft, rectangle));

end
