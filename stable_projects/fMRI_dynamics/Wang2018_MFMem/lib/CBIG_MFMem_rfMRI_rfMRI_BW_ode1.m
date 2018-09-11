
function dF = CBIG_MFMem_rfMRI_rfMRI_BW_ode1(y,F,Nnodes)

%-----------------------------------------------------------------------------
% dF = CBIG_MFMem_rfMRI_rfMRI_BW_ode1(y,F,Nnodes)
%
% Function for hemodynamic model diffiential equation 1st order 
%
% Input:
%     - y:        current neural activiy 
%     - F:        current hemodynamic state      
%     - Nnodes:   number of brain regions
%
% Output:
%     - dF:       change in hemodynamic state
%
% Reference: 
%     (Friston 2000), Nonlinear responses in fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%----------------------------------------------------------------------------


%parameters
beta = 1/0.65; %per s
gamma = 1/0.41; %per s
tau = 0.98; %s
alpha = 0.33;
p = 0.34;

dF = zeros(Nnodes,4);  %F-> Nodes x 4states[dz,df,dv,dq]

dF(:,1) = y - beta*F(:,1) - gamma*(F(:,2)-1); % dz: signal
dF(:,2) = F(:,1);                             % df: flow 
dF(:,3) = 1/tau*(F(:,2)-F(:,3).^(1/alpha));   % dv: volume
dF(:,4) = 1/tau*(F(:,2)/p.*(1-(1-p).^(1./F(:,2)))-F(:,4)./F(:,3).*F(:,3).^(1/alpha)); %dq: deoxyHb


end


