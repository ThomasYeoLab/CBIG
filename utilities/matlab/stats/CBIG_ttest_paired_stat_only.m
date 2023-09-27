function t = CBIG_ttest_paired_stat_only(x)

% Calculate one sample t test value.
%
%     t = CBIG_ttest_paired_stat_only(x)
%     Input:
%         x: N x M matrix, N is number of samples, M is number of trails.
%     Output:
%         t: 1 x M vector, M is number of trails
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if (find(size(x) == 1) == 1)
    x = x';
end

t = mean(x); 
n = size(x,1); 
s = std(x);
t = t./(s./sqrt(n)); 
