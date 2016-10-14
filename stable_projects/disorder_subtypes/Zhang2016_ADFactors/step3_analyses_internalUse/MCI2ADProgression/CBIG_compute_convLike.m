function convLike = CBIG_compute_convLike(b, subtype_compo, gender, avg_stats)

% convLike = CBIG_compute_convLike(b, subtype_compo, gender, avg_stats)
% 1 for male, 2 for female

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

avg_stats(2) = gender;

X = [1 subtype_compo(1:2) avg_stats];
OR = exp(X*b); % log(µ/(1-µ)) = Xb

convLike = OR/(1+OR);
