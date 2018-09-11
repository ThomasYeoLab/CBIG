function dy = CBIG_MFMem_rfMRI_mfm_ode1(y,parameter,SC)

%-----------------------------------------------------------------------------
% dy = CBIG_MFMem_rfMRI_MFMem_ode1(y,parameter,SC)
%
% Function for dynamical mean field model diffiential equation 1st order 
%
% Input:
%     - y:        current neural state
%     - paremter: model parameter vector {p x 1}, in order [w;I0;G], w: self-connection, I0: background, G: global scaling of SC   
%     - SC:       structural connectivity matrix 
%
% Output:
%     - dy:       change in neural state
%
% Reference: 
%     (Deco 2013), Resting-state functional connectivity emerges from structurally and dynamically shaped slow linear fluctuations.
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%----------------------------------------------------------------------------


%parameters for inputs and couplings
J = 0.2609; %nA
w = parameter(1:size(SC,1));
G = parameter(end-1);
I0 = parameter(1+size(SC,1):2*size(SC,1));

%parameters for firing rate
a = 270; %pC
b = 108; %kHz
d = 0.154;   %ms

%parameters for synaptic activity/currents
tau_s = 0.1;  %s 
gamma_s = 0.641;% 


%% total input x

x = J*w.*y+J*G*SC*y+I0;

%% firing rate

c = (1-exp(-d*(a*x-b)));
if c == 0
   disp('error, check firing rate function')
else 
    H = (a*x-b)./c;
end


%% synaptic activity / currents

dy = -1/tau_s*y + gamma_s*(1-y).*H;


end
