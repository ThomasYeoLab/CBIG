function p = CBIG_corrected_resampled_ttest_kfoldCV(accuracy_vec, threshold)

% p = CBIG_corrected_resampled_ttest_kfoldCV(accuracy_vec, threshold)
% 
% This function tests whether an accuracy vector (obtained from repeated K-fold cross-validation) 
% is statistically different from a threshold. The accuracy vector can also
% represent the difference between two classifiers. The key point of this
% function is that it tries to correct for the dependencies between the
% folds, because the sample variance across folds underestimates the true
% variance.
% 
% Input:
%     - accuracy_vec = R x K vector, where R is the number of repetitions and K is the number of folds. 
%
%     - threshold    = 1 x 1 scalar. Null hypothesis: mean(accuracy_vec) 

% Output:
%     - p            = 1 x 1 scalar. 2-sided p value
%
%
% Example:
%     - p = CBIG_corrected_resampled_ttest_kfoldCV([0.5, 0.6, 0.9, 0.4, 0.6], 0.5)
%       This tests whether the accuracies of the 1-time 5-fold CV is statistically
%       different from 0.5.
%
%     - 5-fold CV repeated 2 times: acc = [0.5, 0.6, 0.9, 0.4, 0.6; 0.4, 0.3, 0.2, 0.1, 0.6];
%       p = CBIG_TRBPC_corrected_resampled_ttest(acc, 0);
%       This tests whether the accuracies of the 5-fold CV repeated 2 times is statistically
%       different from 0.
%
%     - Comparing if accuracy numbers given by two methods are significantly
%       different. Assuming for each method we perform 20-fold cross-validation
%       10 times, we will have a accuracy matrix of 10*20 for each model (acc1 and acc2).
%       p = CBIG_TRBPC_corrected_resampled_ttest(acc1-acc2, 0);
%       This tests whether the accuracy of two methods are significantly different.
%
% Reference: 
% 1) Nadeau C, Bengio Y. Inference for the generalization error. NIPS, 2000.
%
% 2) Bouckaert RR, Frank E. Evaluating the replicability of significance tests 
%    for comparing learning algorithms. PAKDD, 2004.
%
% Written by Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


[~, K] = size(accuracy_vec);
assert(K > 1, 'Number of folds must be greater than 1');

% compute the portion of test subjects / training subjects
portion = (1/K) / (1 - 1/K);
p = CBIG_corrected_resampled_ttest(accuracy_vec(:), portion, threshold);