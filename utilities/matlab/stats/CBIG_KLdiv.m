function div = CBIG_KLdiv(ref_prob, prob)

% Calculate KL divergence.
%
%     div = CBIG_KLdiv(ref_prob, prob)
%     Input: 
%         ref_prob: reference probability, Nx1 vector
%         prob    : probability, Nx1 vector
%     Output:
%         div     : KL divergence
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(ref_prob,2) ~= 1
    error('Input argument ''ref_prob'' should be a column vector');
end
if size(prob,2) ~= 1
    error('Input argument ''prob'' should be a column vector');
end

div = sum(ref_prob .* (log2(ref_prob) - log2(prob)));
