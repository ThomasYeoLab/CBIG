function [P,p,Em,En,EN] = caret_P(c,k,u,df,STAT,R)
% Returns the (un)corrected P value using the expected Euler characteristic
% FORMAT [P p Em En EN] = caret_P(c,k,u,df,STAT,R)
%
% c     - cluster number 
% k     - extent {RESELS}
% u     - height-threshold
% df    - [df{interest} df{error}]
% STAT  - Nature of Statisical field
%		'Z' - Gaussian feild
%		'T' - T - feild
%		'X' - Chi squared feild
%		'F' - F - feild
%		'P' - Posterior probability
% R     - RESEL Count {defining search volume}
%
% P     - corrected   P value  - P(n > kmax}
% p     - uncorrected P value  - P(n > k}
% Em    - expected total number of maxima {m}
% En    - expected total number of resels per cluster {n}
% EN    - expected total number of voxels {N}
%
% caret_P returns the probability of c or more clusters with more than
% k surface area (RESELS) thresholded at u.  All p values
% can be considered special cases:
%
% caret_P(1,0,Z,df,STAT,1) = uncorrected p value
% caret_P(1,0,Z,df,STAT,R) = corrected p value {based on height Z)
% caret_P(1,k,u,df,STAT,R) = corrected p value {based on extent k at u)
% caret_P(c,k,u,df,STAT,R) = corrected p value {based on number c at k and u)
% caret_P(c,0,u,df,STAT,R) = omnibus   p value {based on number c at u)
%
% Ref: Hasofer AM (1978) Upcrossings of random fields
% Suppl Adv Appl Prob 10:14-21
% Ref: Friston et al (1993) Comparing functional images: Assessing
% the spatial extent of activation foci
% Ref: Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
% This function is heavily based on spm_P 2.12 Karl Friston 02/10/16
%___________________________________________________________________________
% caret_P.m	1.0 Joern Diedrichsen 04/11/18

%---------------------------------------------------------------------------
% get EC densities
D       = max(find(R))-1;
R       = R(1:D+1);
G       = sqrt(pi)./gamma(([1:D+1])/2);
EC      = spmj_ECdensity(STAT,u,df);
EC      = EC([1:D+1])';


%---------------------------------------------------------------------------
% Get expected values 
EM      = (R.*EC);			% <maxima> over D dimensions
Em      = sum(EM);			% <maxima>
EN      = EC(1)*R(D+1);     % Number of activated resels = EC * activated area
En      = EN/Em;			% Expected size of clusters 

%---------------------------------------------------------------------------
% Now assume that m ~ Poisson and n ~ Exp 
% get p(n>k)
if     ~k | ~D
	p    = 1;
else
	p    = exp(-k/En);
end

%---------------------------------------------------------------------------
% Poisson clumping heuristic {for multiple clusters}
P       = 1 - spm_Pcdf(c - 1,(Em + eps)*p);


% set P and p = [] for non-implemented cases
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if k > 0 & (STAT == 'X' | STAT == 'F')
	P    = []; p = [];
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
