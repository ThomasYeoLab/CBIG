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
%     - One time 5-fold CV: acc = [0.5, 0.6, 0.9, 0.4, 0.6];
%       p = CBIG_TRBPC_corrected_resampled_ttest(acc, 0.25, 0.5);
%       This tests whether the accuracies of the 1-time 5-fold CV is statistically
%       different from 0.5. Note the portion is 1/4 = 0.25 not 1/5.
%
%     - 5-fold CV repeated 2 times: acc = [0.5, 0.6, 0.9, 0.4, 0.6; 0.4, 0.3, 0.2, 0.1, 0.6];
%       p = CBIG_TRBPC_corrected_resampled_ttest(acc(:), 0.25, 0);
%       This tests whether the accuracies of the 5-fold CV repeated 2 times is statistically
%       different from 0. Prior to computing the p-value, we need to
%       vectorize the 2*5 accuracy matrix to a 10*1 vector by using acc(:);
%
%     - train / test split repeated 100 times, each time 70% subjects are used for
%       traning and 30% are used for testing. acc is a 100*1 vector.
%       p = CBIG_TRBPC_corrected_resampled_ttest(acc, 3/7 , 0.6);
%       This tests whether the accuracy is significantly different than 0.6.
%
%     - Comparing if accuracy numbers given by two methods are significantly
%       different. Assuming for each method we perform 20-fold cross-validation
%       10 times, we will have a accuracy matrix of 10*20 for each model (acc1 and acc2).
%       p = CBIG_TRBPC_corrected_resampled_ttest(acc1(:)-acc2(:), 1/19, 0);
%       This tests whether the accuracy of two methods are significantly different.
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
corrected_variance = (1/K + portion) * var(accuracy_vec);

% tstat
mu = mean(accuracy_vec);
tval = (mu-threshold) / sqrt(corrected_variance);

% 2-tail p value (degree of freedom is K - 1)
p = 2 * tcdf(-abs(tval), K-1);