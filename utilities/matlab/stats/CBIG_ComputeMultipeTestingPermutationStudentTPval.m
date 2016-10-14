function val = CBIG_ComputeMultipeTestingPermutationStudentTPval(vec1, vec2, one_sample, N, alpha)

% One or two sample t-test via permutation.
%   
%   val = CBIG_ComputeMultipeTestingPermutationStudentTPval(vec1, vec2, one_sample, N, alpha)
%   Input:
%       vec1        : M1 x T where M is number of subjects and T is the number of tests
%       vec2        : M2 x T where M is number of subjects and T is the number of tests
%       one_sample  : 1-one sample (ttest); 0-two sample (ttest2)
%       N           : number of resamples
%       alpha       : significance lever
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(~exist('one_sample', 'var'))
   one_sample = 1; 
end

if(~exist('N', 'var'))
   N = 1000; 
end

if(~exist('alpha', 'var'))
   alpha = 0.05; 
end

if(one_sample)
    diff_size = size(vec1) - size(vec2);
    if(sum(abs(diff_size)) > 0)
        error('vec1 and vec2 not the same size!');
    end
else
    if(size(vec1, 2) ~= size(vec2, 2))
        error('Number of tests different!');
    end
end

rand('state',sum(100*clock));

combinedVec = [vec1; vec2];
M1 = size(vec1, 1);
M2 = size(vec2, 1);
TotalSubjects = size(combinedVec, 1);
val = zeros(N, 1);
for i = 1:(N-1)
    
    newVec = combinedVec(randperm(TotalSubjects), :);
    if(one_sample)
        newVecDiff = newVec(1:end/2, :) - newVec(end/2+1:end, :);
        [~, b] = ttest(newVecDiff);
    else
        [~, b] = ttest2(newVec(1:M1, :), newVec(M1+1:end, :));
    end
    val(i) = min(b);
end

w = sort(val, 'ascend');
val = w(round(alpha*N));
