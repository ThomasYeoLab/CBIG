function p = CBIG_corrected_resampled_ttest(accuracy_vec, portion, threshold)

% p = CBIG_corrected_resampled_ttest(accuracy_vec, portion, threshold)
% 
% This function tests whether an accuracy vector (obtained from cross-validation) 
% is statistically different from a threshold. The accuracy vector can also
% represent the difference between two classifiers. The key point of this
% function is that it tries to correct for the dependencies between the
% folds, because the sample variance across folds underestimates the true
% variance.
% 
% Input:
%     - accuracy_vec = K x 1 vector, K is the number of folds. 
%
%     - portion      = 1 x 1 scalar. # test subjects/# training subjects
%                                    for each cross-validation folds. E.g. portion for
%                                    5-fold CV equals 1/4 = 0.25
%     - threshold    = 1 x 1 scalar. Null hypothesis: mean(accuracy_vec) 

% Output:
%     - p            = 1 x 1 scalar. 2-sided p value
%
%
% Example:
%     - p = CBIG_TRBPC_corrected_resampled_ttest([0.5, 0.6, 0.9, 0.4, 0.6], 0.25, 0.5)
%       This tests whether the accuracies of the 1-time 5-fold CV is statistically
%       different from 0.5.
%
% Reference: 
% 1) Nadeau C, Bengio Y. Inference for the generalization error. NIPS, 2000.
%
% 2) Bouckaert RR, Frank E. Evaluating the replicability of significance tests 
%    for comparing learning algorithms. PAKDD, 2004.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


K = length(accuracy_vec);
assert(K > 1, 'Number of folds must be greater than 1');

% corrected variance
% n2 = # samples in test set, n1 = 1 / K (for K-fold CV)
% n1 = # samples in training set, n2 = (1 - 1/K) (for K-fold CV)
corrected_variance = (1/K + portion) * var(accuracy_vec);

% tstat
mu = mean(accuracy_vec);
tval = (mu-threshold) / sqrt(corrected_variance);

% 2-tail p value (degree of freedom is K*R - 1)
p = 2 * tcdf(-abs(tval), K-1);
