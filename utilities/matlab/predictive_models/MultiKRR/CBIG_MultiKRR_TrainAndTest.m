function [pred_acc,corr_acc,pred_stats, y_predicted] = CBIG_MultiKRR_TrainAndTest(FSM_train,...
    FSM_pred,y_train_resid,y_pred_resid,y_pred_orig,lambda_vect,...
    bin_flag,with_bias,threshold, acc_metric, train_test_flag)

% [pred_acc,corr_acc,pred_stats, y_predicted] = CBIG_MultiKRR_TrainAndTest(FSM_train,...
%   FSM_pred,y_train_resid,y_pred_resid,y_pred_orig,num_kernel,lambda_vect,...
%   bin_flag,with_bias,threshold, acc_metric)
%
% This function performs the training and prediction for one behavioral
% measure using multi KRR model. Note that for training, the function will
% use the selected hyperparameter and output an accuracy value based on
% what the user specified in `acc_metric`. For testing, the function will
% output all possible acc_metric values given the hyperparameters into a
% variable called `pred_stats`.
% 
% Inputs:
%   - FSM_train
%     A #train_sub by #train_sub by #kernel matrix. 
%
%   - FSM_pred
%     A #test_sub by #train_sub by #kernel matrix.
%
%   - y_train_resid
%     A #train_sub by 1 vector. This contains the regressed behavioral
%     measure in the training set.
%     
%   - y_pred_resid
%     A #test_sub by 1 vector. This contains the regressed behavioral
%     measure in the test set.
%
%   - y_pred_orig
%     A #test_sub by 1 vector. This contains the original pre-regressed
%     behavioral measure. This is only necessary if the behavioral measure
%     is binary. If not, the user can pass in ''.
%
%   - lambda_vect
%     A 1 by #kernel vector. lambda_vect(1,i) stores the hyperparamter used to 
%     regularize the ith kernel (ie. FSM(:,:,i)). 
%
%   - bin_flag
%     A scalar. If bin_flag is 1, this means the behavioral measure is
%     binary. If bin_flag is 0, this means the behavioral measure is
%     continuous.
%
%   - with_bias
%     A binary scalar (choose between 0 or 1).
%     - with_bias = 0 means the algorithm is to minimize 
%     (y - K*alpha)^2 + (regularization of alpha);
%     - with_bias = 1 means the algorithm is to minimize
%     (y - K*alpha - beta)^2 + (regularization of alpha), where beta is a
%     constant bias for every subject, estimated from the data.
%
%   - threshold
%     A scalar, it is a threshold required to determine the prediction (1
%     or 0) for binary target variables (e.g. sex). 
%
%   - acc_metric
%     A string stating which accuracy metric to be used for optimising
%     hyperparameters. Currently, the following are accepted:
%     Choose from:
%       'corr'              - Pearson's correlation;
%       'COD'               - Coefficient of determination. Defined as
%                             1-||y_pred-y_test||^2/||mean(y_test)-y_test||^2,
%                             where y_pred is the prediction of the test data, 
%                             y_test is the groud truth of the test data, 
%                             and mean(y_test) is the mean of test data
%       'predictive_COD'    - Predictive coefficient of determination. Defined as
%                             1-||y_pred-y_test||^2/||mean(y_train)-y_test||^2,
%                             where y_pred is the prediction of the test data, 
%                             y_test is the groud truth of the test data, 
%                             and mean(y_train) is the mean of training data
%       'MAE'               - mean absolute error
%       'MAE_norm'          - mean absolute error divided by the standard
%                             derivation of the target variable of the training set
%       'MSE'               - mean squared error
%       'MSE_norm'          - mean squared error divided by the variance
%                             of the target variable of the training set
%
%   - train_test_flag
%     A scalar stating if the user wishes to perform training or testing.
%     This is necessary because during test, this function will tabulate all the
%     possible acc_metric results (ie 'corr','COD', ....'MSE_norm') based
%     on the optimised hyperparameter trained on acc_metric.0 means train
%     and 1 means test.
% 
% Outputs:
%   - pred_acc
%     A scalar. Stores the prediction accuracy (negative loss) after training 
%     and testing.The metric used for calculating the pred_acc is
%     determined by acc_metric.
%
%   - corr_acc
%     A scalar. Stores the prediction accuracy after training 
%     and testing.The metric used for calculating the pred_acc is
%     correlation
%
%   - pred_stats
%     A data structure. There are 2 fields in the data structure, 'value' and
%     'description'. 'description' describes which accuracy metric the
%     'value' belongs to. 'value' will be a scalar
%     storing the values of the optimal test prediction accuracies
%     corresponding to the specified acuracy metric. Within each field,
%     there will be a total number of 7 entries, each entry corresponding
%     to the possible acc_metrics.
%
%   - y_predicted
%     A #test_sub by 1 vector. Stores the predicted values after trainin
%     and testing using the multi KRR model. 
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
% Author: Yan Rui Tan and Ru(by) Kong

%% Check input vectors
if size(y_train_resid,2) > 1
    error('"y_train_resid" should be a column vector')
end

if size(y_pred_resid,2) > 1
    error('"y_pred_resid" should be a column vector')
end

if size(y_pred_orig,2) > 1
    error('"y_pred_orig" should be a column vector')
end

if size(lambda_vect,1) > 1
    error('"lambda_vect" should be a row vector')
end

if ischar(train_test_flag)
    train_test_flag = str2num(train_test_flag);
end

%% Check nan index and remove from subjects
idx_nan_train = ~isnan(y_train_resid);
idx_nan_pred = ~isnan(y_pred_resid);
FSM_train = FSM_train(idx_nan_train, idx_nan_train,:);
FSM_pred = FSM_pred(idx_nan_pred, idx_nan_train,:);
y_train_resid = y_train_resid(idx_nan_train, 1);
y_pred_resid = y_pred_resid(idx_nan_pred, 1);
y_pred_orig = y_pred_orig(idx_nan_pred, 1);

%% Determine number of training and test subjects and number of kernels used
num_kernel = size(FSM_train,3);
num_train_subj = size(FSM_train,1);
num_test_subj = size(FSM_pred,1);

%% Perform training and prediction
y_tmp = repmat(y_train_resid,num_kernel,1);
K_tmp = repmat(reshape(FSM_train,num_train_subj,num_train_subj*num_kernel),...
    num_kernel,1);
kr_eye_tmp = eye(num_kernel*num_train_subj);
kr_lambda_tmp = reshape(repmat(lambda_vect,num_train_subj,1),...
    num_train_subj*num_kernel,1);
kr_tmp = bsxfun(@times, kr_eye_tmp, kr_lambda_tmp);
K_tmp = K_tmp + kr_tmp;

if with_bias == 0
    %Training
    alpha = K_tmp \ y_tmp;
    %Prediction
    y_predicted = reshape(FSM_pred,num_test_subj,num_train_subj*num_kernel)*alpha;
else
    %Training
    inv_K_tmp = inv(K_tmp);
    X = ones(num_train_subj*num_kernel,1);
    beta = (X' * (inv_K_tmp * X)) \ X' * (inv_K_tmp * y_tmp);
    alpha = inv_K_tmp * (y_tmp - X * beta);
    %Prediction
    y_predicted = reshape(FSM_pred,num_test_subj,...
        num_train_subj*num_kernel)*alpha + ones(num_test_subj,1) .* beta;
end

%% Find prediction accuracies
if bin_flag ==1
    TP = length(find((y_predicted > threshold) & (y_pred_orig==1) == 1));
    TN = length(find((y_predicted < threshold) & (y_pred_orig==0) == 1));
    pred_acc = (TP + TN) / length(y_pred_orig);
    corr_acc = nan;
else
    if train_test_flag == 0
        [pred_stats,loss] = CBIG_compute_prediction_acc_and_loss(y_predicted,y_pred_resid,...
            acc_metric,y_train_resid);
        pred_acc = -loss;
        corr_acc = nan;
       
    else
        all_acc_metric = {'corr', 'COD', 'predictive_COD','MAE', 'MAE_norm',...
            'MSE','MSE_norm'};
        for jj = 1:length(all_acc_metric)
            curr_acc_metric = all_acc_metric{jj};
            [pred_stats,loss] = CBIG_compute_prediction_acc_and_loss(y_predicted,y_pred_resid,...
            curr_acc_metric,y_train_resid);
            if jj == 1
                corr_acc = pred_stats;
            end
            opt_stats(jj).description= curr_acc_metric;
            opt_stats(jj).value = pred_stats;
        end
        pred_stats = opt_stats;
        pred_acc = nan;
    end
            
end

end
