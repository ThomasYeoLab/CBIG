function CBIG_RL2017_Figure4_toy_ex_HMM

% CBIG_RL2017_Figure4_toy_ex_HMM
% This code generates a 2-state HMM model and the corresponding autoregressive (AR) and phase 
% randomized (PR) surrogates. It was used to generate Figure 4 in Liegeois et al. (NeuroImage 2017).
%
%
% Written by Raphael Liegeois and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Contact: Raphael Liegeois (Raphael.Liegeois@gmail.com)


% Generating original time course
N = 1201; %to be consistent with HCP fMRI time series

% Defining two states
Corr1 = [1 -0.2; -0.2 1]; R1 = chol(Corr1);
Corr2 = [1 0.9; 0.9 1];   R2 = chol(Corr2);

% Tranistion probability matrix
Transition = [0.995 0.005;0.005 0.995];

State = 2; % Initial state, without loss of generality

x = zeros(N,2);
Transition_prob = rand(1,N); %Used to determine when to change state

for i=1:N
    if State==1
        x(i,:)=randn(1,2)*R1;
        if Transition_prob(i)>Transition(1,1)
            State=2; %change state
        end
    else
        x(i,:)=randn(1,2)*R2;
        if Transition_prob(i)>Transition(2,2)
            State=1; %change state
        end
    end
end

% Compute AR and PR surrogates
x_AR = CBIG_RL2017_get_AR_surrogate(x,1,1);
x_PR = CBIG_RL2017_get_PR_surrogate(x,1);

% Compute sliding window correlation
w    = 83; %window width, same as for the real fMRI data analysis

for i = 1:N-w
    t_orig         = corr(x(i:i+w,:));
    x_corr(1,i)    = t_orig(1,2); % consider (windowed) correlation between two variables
    t_ar           = corr(x_AR(i:i+w,:));
    x_corr_ar(1,i) = t_ar(1,2);
    t_pr           = corr(x_PR(i:i+w,:));
    x_corr_pr(1,i) =t_pr(1,2);
end

% Plot Results

subplot(1,3,1)
hold on
plot(x(:,1),'b')
plot(x(:,2),'g')
plot([ceil(w/2):ceil(w/2)+length(x_corr)-1],x_corr,'r','linewidth',3) %centering the sliding window correlation
title('Original 2-state HMM','fontsize',20)
axis([0 N-w+1 -5 5]) 

subplot(1,3,2)
hold on
plot(x_AR(:,1),'b')
plot(x_AR(:,2),'g')
plot([ceil(w/2):ceil(w/2)+length(x_corr_ar)-1],x_corr_ar,'r','linewidth',3)
title('AR surrogate','fontsize',20)
axis([0 N-w+1 -5 5])

subplot(1,3,3)
hold on
plot(x_PR(:,1),'b')
plot(x_PR(:,2),'g')
plot([ceil(w/2):ceil(w/2)+length(x_corr_pr)-1],x_corr_pr,'r','linewidth',3)
title('PR surrogate','fontsize',20)
axis([0 N-w+1 -5 5])

