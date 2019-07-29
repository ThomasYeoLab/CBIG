function [obs_wc, unobs_wc] = CBIG_SplitWCIntoObsUnObs(wc, frac_obs)

% Split words into observed and unobserved
% FORMAT [obs_wc, unobs_wc] = SplitWCIntoObsUnObs(wc, frac_obs)
%
% wc       = 1 x V vector, where wc(i) = # times ith vocab word appears
% frac_obs = fraction of words that are deemed observable
%
% obs_wc   = observed words
% unobs_wc = unobservved words 
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(wc,1) ~= 1
    error('Input argument ''wc'' should be a row vector');
end

V = size(wc, 2);
[~, I] = sort(rand(1, V), 2);

wc_permute = wc(I);
wc_permute_cumsum = cumsum(wc);
wc_threshold = frac_obs*sum(wc);

threshold = find(wc_permute_cumsum > wc_threshold, 1);


obs_wc = zeros(1, V);
obs_wc(I(1:threshold)) = wc(I(1:threshold));

unobs_wc = zeros(1, V);
unobs_wc(I(threshold+1:end)) = wc(I(threshold+1:end));


