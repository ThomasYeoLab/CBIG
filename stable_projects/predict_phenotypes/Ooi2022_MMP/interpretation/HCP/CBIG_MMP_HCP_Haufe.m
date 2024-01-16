function CBIG_MMP_HCP_Haufe(input_dir, results_dir, feature)

% function CBIG_MMP_HCP_Haufe(input_dir, results_dir, feature)
%
% This function calculated the Haufe-inverted feature importance for each model. This function
% is specifically for KRR models and not meant to invert the features from LRR or Elasticnet.
%
% Input:
% - input_dir
% The directory in which the brain imaging features are saved.
%
% - results_dir
% The directory in which the regression results are results are saved.
%
% - feature
% The outstem of the model to invert (e.g. features_rs).
%
% Output:
% - cov_mat_mean
% A mat file is saved with a matrix of #seeds x #features x #behaviors. Each entry of the matrix represents the
% feature importance averaged over the outerfolds for each seed.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% process for each feature
save_dir = fullfile(results_dir, 'interpretation', feature);
if ~exist(save_dir)
    mkdir(save_dir);
end
input_mat = fullfile(input_dir, strcat(feature, '.mat'));
model = strcat('KRR_', feature);
total_seeds = 60;

% load seed
if exist(fullfile(save_dir, 'cov_mat_mean.mat'))
    fprintf('cov_mat_mean exists. File will not be generated.\n')
else
    clear cov_mat_mean
    if exist(fullfile(save_dir, 'cov_mat_tmp.mat'))
        load(fullfile(save_dir, 'cov_mat_tmp.mat'))
        fprintf('Continuing from previous run, continue from seed %i \n', last_s)
    else
        last_s = 1;
    end
    for s = last_s:total_seeds
        fprintf('Calculating for %s, seed %i / %i \n', feature, s, total_seeds)
        seed_name = strcat('seed_', num2str(s));
        model_dir = fullfile(results_dir, model, seed_name, 'results');
        
        % load sub_fold
        load(fullfile(model_dir, 'no_relative_10_fold_sub_list.mat'));
        % load results
        results = load(fullfile(model_dir, strcat('final_result_', model, '.mat')));
        % load features
        feat = load(input_mat);
        
        % pre allocate space for cov_mat
        cov_mat = zeros(size(sub_fold,1), size(feat.(feature),1), ...
            size(results.y_pred_train{1},2));
        
        % calculate feature importance for each fold
        for i = 1:size(sub_fold,1)
            % find train fold idx
            train = ~sub_fold(i).fold_index;
            fold_name = strcat('fold_', num2str(i));
            % load features and normalize
            feat_train = feat.(feature)(:,train);
            feat_train_norm = (feat_train - mean(feat_train,1)) ./ std(feat_train, [], 1);
            % load predictions
            y_pred = results.y_pred_train{i};
            
            % compute covariance
            for b = 1:size(y_pred,2)
                cov_mat(i,:,b) = bsxfun(@minus,feat_train_norm,mean(feat_train_norm,2)) * ...
                    bsxfun(@minus,y_pred(:,b),mean(y_pred(:,b))) / (size(feat_train_norm,2));
            end
        end
        
        % take average over outerfolds
        cov_mat_mean(s,:,:) = mean(cov_mat,1);
        save(fullfile(save_dir, 'cov_mat_tmp.mat'), 'cov_mat_mean', 'last_s', '-v7.3');
        last_s = s;
    end
    % save
    delete(fullfile(save_dir, 'cov_mat_tmp.mat'))
    save(fullfile(save_dir, 'cov_mat_mean.mat'), 'cov_mat_mean', '-v7.3');
end

end