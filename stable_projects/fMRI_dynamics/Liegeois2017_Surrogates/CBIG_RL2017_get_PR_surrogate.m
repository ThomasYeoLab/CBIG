function surr = CBIG_RL2017_get_PR_surrogate(TC , n_surr)

% surr = CBIG_RL2017_get_PR_surrogate(TC , n_surr)
% This function generates phase randomized (PR) surrogates. The two
% mandatory inputs are:

% 'TC'                 original time course of size (num_timepoints,num_variables)
%                      In this case num_variables is the numbers of ROIs.
% 'n_surr'             desired number of surrogate samples to be generated

% Output
% 'surr'               matrix of size (num_timepoints,num_variables,n_surr)
%                      containing surrogate datasets

% Example
% x      = randn(101,4);
% x(:,4) = sum(x,2); % Create correlation in the data
% r1     = corrcoef(x) ;
% surr   = CBIG_RL2017_get_PR_surrogate(x,2);
% r2     = corrcoef(surr(:,:,2)); % Check that the correlation is preserved
% norm(r1-r2); 

% This code is based on previous work by Carlos Gias.
% Written by Raphael Liegeois and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Contact: Raphael Liegeois (Raphael.Liegeois@gmail.com)

%%% References: 
%               * Prichard, D., Theiler, J. Generating Surrogate Data for Time Series
%                 with Several Simultaneously Measured Variables (1994), Physical Review Letters, Vol 73, Number 7
%               * For review: R. Liegeois et al. (2017). Interpreting Temporal
%                 Fluctuations In Resting-State Functional Connectivity MRI, NeuroImage.

%% Identify parameters
T  = size(TC,1);
k  = size(TC,2);

if rem(T,2)==0 %make the number of samples odd
    T  = T-1;
    TC = TC(1:T,:);
end

len_ser = (T-1)/2;
interv1 = 2:len_ser+1; 
interv2 = len_ser+2:T;

% Fourier transform of the original dataset
fft_TC = fft(TC);

%% Create the surrogate data
surr = zeros(T, k, n_surr);
for s = 1:n_surr
   ph_rnd = rand([len_ser 1]);
   
   % Create the random phases for all the time series
   ph_interv1 = repmat(exp(2*pi*1i*ph_rnd),1,k);
   ph_interv2 = conj(flipud(ph_interv1));
   
   % Randomize all the time series simultaneously
   fft_TC_surr = fft_TC;
   fft_TC_surr(interv1,:) = fft_TC(interv1,:).*ph_interv1;
   fft_TC_surr(interv2,:) = fft_TC(interv2,:).*ph_interv2;
   
   % Inverse transform
   surr(:,:,s)= real(ifft(fft_TC_surr));
end
