function CBIG_RL2017_Figure7_toy_ex_bi_vs_multi

% CBIG_RL2017_Figure7_toy_ex_bi_vs_multi
% 
% This code illustrates the difference between multivariate and bivariate
% approaches to generate surrogate time series. The multivariate approach 
% consists in identifying a single multivariate AR model of the original time 
% series whereas the bivariate AR approach consists in identifying a bivariate
% AR model for each pair of variables. We here illustrate that
% the bivariate approach induces biases in the identification of AR parameters. 
% This code was used to generate Figure 7 in Liegeois et al. (NeuroImage 2017).
%
%
% Written by Raphael Liegeois and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Contact: Raphael Liegeois (Raphael.Liegeois@gmail.com)

% Defining Variables
A = [.3 0.4 -0.5; 0.4 .3 0; -0.5 0 .3]; %True AR parameters
N = 100000; %Number of time points
x = zeros(3,N);

% Checking that the AR model is stable (absolute eigenvalues must be < 1)
e_val=eig(A);

if max(abs(e_val))<1
    disp('Sanity check: AR model is stable')
else
    error('AR model is not stable and the identification procedure might not be meaningful')
end

% Generate Time courses from this AR model
e      = randn(3,N);  %Input noise
x(:,1) = e(:,1); %Initialization

for i=2:N
    x(:,i)=A*x(:,i-1)+e(:,i);
end

% Compute back the AR parameters in the multivariate case

[Y,B_multi,Z,E]=CBIG_RL2017_ar_mls(x',1);

disp('True AR(1) model')
A
disp('Identified multivariate AR model (good match with ground truth):')
B_multi(:,2:4)
disp('Three bivariate AR models:')
% First pairwise identification
x12=x([1 2],:);
[Y,B_12,Z,E]=CBIG_RL2017_ar_mls(x12',1);
B_12(:,2:3)

%second permutation
x13=x([1 3],:);
[Y,B_13,Z,E]=CBIG_RL2017_ar_mls(x13',1);
B_13(:,2:3)

%third permutation
x23=x([2 3],:);
[Y,B_23,Z,E]=CBIG_RL2017_ar_mls(x23',1);
B_23(:,2:3)

disp('Note that there is no link between the 2nd and 3rd toy brain regions in the ground truth A matrix. However, in the bivariate estimation built from these two variables (x2 and x3), the link between the two brain regions is nonzero (~0.25).')

