function [FC_simR, CC_check] = CBIG_MFMem_rfMRI_nsolver_eul_sto_4p(parameter,prior,SC,y_FC,FC_mask,Nstate,Tepochlong,TBOLD,ifRepara)

%-----------------------------------------------------------------------------
% FC_simR, CC_check] = CBIG_mfm_rfMRI_nsolver_eul_sto(parameter,prior,SC,y_FC,FC_mask,Nstate,Tepochlong,TBOLD,ifRepara)
%
% Function to 
%  (a)solve diffitial equation of dynamic mean field and hemodynamic model using stochastic Euler method 
%  (b)caculate simulated functional connectivity (FC) from simulated BOLD
%  (c)caculate correlation between emprical FC and simulated FC 
%
% Input:
%     - SC:        structural connectivity matrix
%     - y_FC:      functional connectivity vector (after mask)
%     - FC_mask:   used to select uppertragular elements
%     - Nstate:    noise randon seed
%     - Tepochlong:simulation long in [min], exclusive 2min pre-simulation(casted)
%     - TBOLD:     BOLD-signal time resolution
%     - ifRepara:  re-parameter turn on=1/off=0
%
% Output:
%     - FC_simR:  simulated FC, only entries above main diagonal, in vector form
%     - CC_check: cross correlation of 2 FCs 
%
% Reference:
%     [1](Deco 2013), Resting-state functional connectivity emerges from structurally and dynamically shaped slow linear fluctuations.
%     [2](Friston 2000), Nonlinear responses in fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%----------------------------------------------------------------------------

%caculate first [BOLD,yT,fT,qT,vT,zT,Time], then caculate FC_stim


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%(a) solve diffitial equation of dynamic mean field and hemodynamic model using stochastic Euler method 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if size(parameter, 2) > 1
    error('Input argument ''parameter'' should be a column vector');
end

if ifRepara == 1
    parameter = exp(parameter).*prior;
end


%-----------------------------------------------------------------------
%initial system
%-----------------------------------------------------------------------
%simulation time
kstart = 0;  %s
Tpre = 60*2; %s
kend = Tpre+60*Tepochlong; %s

dt_l = 0.01;    %s  integration time step

dt = dt_l;  %s  time step for neuro
dtt = 0.01; %s, time step for BOLD

%sampling ratio
k_P = kstart:dt:kend;   
k_PP = kstart:dtt:kend;

%initial
Nnodes = size(SC,1);
Nsamples = length(k_P);
Bsamples = length(k_PP);

%for neural activity y0 = 0
yT = zeros(Nnodes,1);

%for hemodynamic activity z0 = 0, f0 = v0 = q0 =1
zT = zeros(Nnodes,Bsamples);
fT = zeros(Nnodes,Bsamples);
fT(:,1) = 1;
vT = zeros(Nnodes,Bsamples);
vT(:,1) = 1;
qT = zeros(Nnodes,Bsamples);
qT(:,1) = 1;

F = [zT(:,1) fT(:,1) vT(:,1) qT(:,1)];
yT(:,1) = 0.001;


%wiener process
w_coef = parameter(end)/sqrt(0.001); %0.001 0.006
w_dt = dt; %s
w_L = length(k_P);
rng(Nstate);
dW = sqrt(w_dt)*randn(Nnodes,w_L+1000);  %plus 1000 warm-up

j = 0;

%--------------------------------------------------------------------------
%solver: Euler
%-------------------------------------------------------------------------- 
tic

%% warm-up


for i = 1:1000
        
        dy = CBIG_MFMem_rfMRI_mfm_ode1_4p(yT,parameter,SC);
        yT = yT + dy*dt + w_coef*dW(:,i);

end


%% main body: caculation 

for i = 1:length(k_P)
        
        dy = CBIG_MFMem_rfMRI_mfm_ode1_4p(yT,parameter,SC);
        yT = yT + dy*dt + w_coef*dW(:,i+1000);
        if mod(i,dtt/dt) == 0
            j = j+1;
            y_neuro(:,j) = yT;
        end
        
end

for i = 2:length(k_PP)

            dF = CBIG_MFMem_rfMRI_rfMRI_BW_ode1(y_neuro(:,i-1),F,Nnodes);
            F = F + dF*dtt;
            zT(:,i) = F(:,1);
            fT(:,i) = F(:,2);
            vT(:,i) = F(:,3);
            qT(:,i) = F(:,4);

end


%Parameter for Balloom-Windkessel model, we updated  this model and its
%parameter according to Stephan et al 2007, NeuroImage 38:387-401 and
%Heinzle et al. 2016 NeuroImage 125:556-570, Parameter are set for the 3T
%and TE=0.0331s
p = 0.34; 
v0 = 0.02;
k1 = 4.3*28.265*3*0.0331*p;
k2 = 0.47*110*0.0331*p;
k3 = 1-0.47;
y_BOLD = 100/p*v0*( k1*(1-qT) + k2*(1-qT./vT) + k3*(1-vT) );    

Time = k_PP;
toc
 


%--------------------------------------------------------------------------
%(b)&(c) compute simulated FC and correlation of 2 FCs
%--------------------------------------------------------------------------
%get the static part
cut_indx = find(Time == Tpre);% after xx s
BOLD_cut = y_BOLD(:,cut_indx:end);
y_neuro_cut = y_neuro(:,cut_indx:end);
Time_cut = Time(:,cut_indx:end);

BOLD_d = CBIG_MFMem_rfMRI_simBOLD_downsampling(BOLD_cut,TBOLD/dtt); %down sample 

FC_sim = corr(BOLD_d');
FC_simR = FC_sim(~FC_mask);


CC_check = corr(atanh(FC_simR),atanh(y_FC));


end

