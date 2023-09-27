function t = CBIG_ttest2_stat_only(x,y)

% Given two distributions, calculate 2 sample t test value.
% 
%     t = CBIG_ttest2_stat_only(x,y)
%     Input:
%         x: N1 x M matrix, N1 is number of samples, M is number of trails.
%         y: N2 x M matrix, N2 is number of samples, M is number of trails.
%     Output:
%         t: 1 x M vector, M is number of trails.
% PS: If x, y is 1 x N or N x 1 vector, the function can handle both cases.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if (find(size(x) == 1) == 1)
    x = x';
end
if (find(size(y) == 1) == 1)
    y = y';
end
if (size(x,2) ~= size(y,2))
    error('ERROR: input matrix x and y should have same width');
end

t = mean(x)-mean(y); 
n1 = size(x,1); 
n2 = size(y,1); 
s = sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2)); 
t = t./(s*sqrt(1/n1+1/n2)); 

if (find(size(x) == 1) == 1)
    t = t';
end

