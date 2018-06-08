function p = CBIG_corrected_resampled_ttest_kfoldCV(accuracy_vec, threshold)

% p = CBIG_corrected_resampled_ttest_kfoldCV(accuracy_vec, threshold)
% 
% This function tests whether an accuracy vector (obtained from K-fold cross-validation) 
% is statistically different from a threshold. The accuracy vector can also
% represent the difference between two classifiers. The key point of this
% function is that it tries to correct for the dependencies between the
% folds, because the sample variance across folds underestimates the true
% variance.
% 
% Input:
%     - accuracy_vec = K x 1 vector, where K is the number of folds. 
%
%     - threshold    = 1 x 1 scalar. Null hypothesis: mean(accuracy_vec) 

% Output:
%     - p            = 1 x 1 scalar. 2-sided p value
%
%
% Example:
%     - p = CBIG_corrected_resampled_ttest_kfoldCV([0.5; 0.6; 0.9; 0.4; 0.6], 0.5)
%       This tests whether the accuracies of the 5-fold CV is statistically
%       different from 0.5.
%
% Reference: 
% 1) Nadeau C, Bengio Y. Inference for the generalization error. NIPS, 2000.
%
% 2) Bouckaert RR, Frank E. Evaluating the replicability of significance tests 
%    for comparing learning algorithms. PAKDD, 2004.
%
% Written by Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% num folds
K = length(accuracy_vec);

% corrected variance
% n2 = # samples in test set, n1 = 1 / K (for K-fold CV)
% n1 = # samples in training set, n2 = (1 - 1/K) (for K-fold CV)
n2_div_n1 = (1/K) / (1 - 1/K);
corrected_variance = (1/K + n2_div_n1) * var(accuracy_vec);

% tstat
mu = mean(accuracy_vec);
tval = (mu-threshold) / sqrt(corrected_variance);

% 2-tail p value (degree of freedom is K - 1)
p = 2 * tcdf(-abs(tval), K-1);