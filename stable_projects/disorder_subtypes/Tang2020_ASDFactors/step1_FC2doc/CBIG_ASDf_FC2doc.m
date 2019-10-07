function [Z, discretized_Z] = CBIG_ASDf_FC2doc(corr_mat, regressors, dx_info, cohort_label, output_dir, output_name)
% [Z, discretized_Z] = CBIG_ASDf_FC2doc(corr_mat, regressors, dx_info, cohort_label, output_dir, output_name)
% 
% Convert FC matrix to "documents" to be used as inputs to polarLDA
% model for model parameters estimation. This function will:
% 1) estimate GLM parameters with only control participants to do regression; 
% 2) compute post-regression control group mean and std for Z-normalization;
% 3) perform Z-normalization with respect to control participants;
% 4) times 10 to Z-score and discretize to integer values so that positive means 
% hyper-connectivity (with respect to controls) and negative means hypo-connectivity; 
% 5) write into documents which can be used as inputs to polarLDA model.
% 
% Input:
%     - corr_mat: 
%           MxMxN FC matrix, where M is the number of ROIs, N is the number of subjects
%     - regressors: 
%           NxV matrix, where N is the number of subjects, V is the num of
%           regressors. The subject's order should be the same as that in corr_mat 
%     - dx_info: 
%           Nx1 vector, where N is the num of subjects. This is the diagnostic information
%           of the subjects. E.g., in ABIDE: 1 represents ASD and 2 represents controls.
%           The subject's order should be the same as that in corr_mat
%     - cohort_label: 
%           Kx1 vector, where K is the number of different cohorts (assuming
%           1st entry of cohort_label represents controls). This function
%           outputs "document" for each cohort.
%           E.g., following ABIDE convention, cohort_label should be [2; 1].
%     - output_dir:
%           Absolute path to the output directory
%     - output_name: 
%           A string to be used as the prefix to the output file name
% Output:
%     - Z:
%           NxP matrix, where P = (Mx(M-1))/2 is the number of unique
%           ROI-ROI pairs. This is the Z-scores after regression and Z-normalization with
%           respect to controls.
%     - discretized_Z:
%           NxP matrix. This is discretized Z-scores and is written into
%           the output documents output_name_dx1.dat and
%           output_name_dx2.dat.
%     - output_name_dx1.dat, output_name_dx2.dat etc.:
%           "Document" of each cohort, saved in the directory specified by output_dir.
%           E.g., output_name_dx1.dat is the "document" of cohort with
%           cohort label 1 (e.g., ASD cohort in the case of ABIDE). This 
%           document can be used as an input to polarLDA model.
%     - output_name_mean_CN.mat:
%           Mean of FC data of control participants, saved in the directory
%           specified by output_dir
%     - output_name_beta_CN.mat:
%           Regression coefficients estimated from control participants,
%           saved in the directory specified by output_dir
%     - output_name_reg_CN_mean_std.mat:
%           Post-regression mean and std of control participants, saved in
%           the directory specified by output_dir
%           
% Example:
%     [Z, discretized_Z] = CBIG_ASDf_FC2doc(corr_mat, regressors, dx_info,
%     [2; 1], '~/step1_FC2doc/output/', 'step1_output')
%     This example computes the discretized Z-scores and writes into two
%     separate documents (for ASD and control cohorts, respectively)
%     in the directory '~/step1_FC2doc/output/' with prefix file name
%     'step1_output', and gives Z-scores and discretized Z-scores as outputs.
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%%Check input variables
if size(dx_info,2) > 1
    error('Input argument ''dx_info'' should be a column vector');
end

if size(cohort_label,2) > 1
    error('Input argument ''cohort_label'' should be a column vector');
end

%% Get number of subjects
no_subjects = size(corr_mat,3);

%% Get unique ROI-ROI correlation values (lower triangle of the FC matrix)
for num = 1:no_subjects
    corr_one_sub = corr_mat(:,:,num);
    lowerTri = tril(ones(size(corr_one_sub)),-1);
    corrArr(num,:) = corr_one_sub(lowerTri ~= 0);
end

%% Regress out the regressors
% Estimate GLM parameters with only control subjects
indCN = dx_info==cohort_label(1);
Y_CN = corrArr(indCN,:);
mean_CN = mean(Y_CN,1);
save([output_dir output_name '_mean_CN.mat'], 'mean_CN');
X_CN = [ones(size(Y_CN,1),1) regressors(indCN,:)];

% add a diagonal matrix with small values to prevent singular matrix problem
b_CN = (X_CN'*X_CN + eye(size(X_CN,2))*1e-6)\(X_CN'*Y_CN);
%b_CN = (X_CN'*X_CN)\(X_CN'*Y_CN);
save([output_dir output_name '_beta_CN.mat'], 'b_CN');

% Regress out regressors for all subjects by computing Y - X*b_CN +
% mean(Y_CN) for each unique ROI-ROI pair
Y = corrArr;
X = [ones(size(Y,1),1) regressors];
corrArr_reg = bsxfun(@plus, Y-X*b_CN, mean_CN);

disp('Regression done.');

%% Z-normalization and discretization
% Compute post-regression control group mean and std for Z-normalization
reg_CN_mean = mean(corrArr_reg(indCN,:),1);
reg_CN_std = std(corrArr_reg(indCN,:),0,1);
save([output_dir output_name '_reg_CN_mean_std.mat'], 'reg_CN_mean', 'reg_CN_std');

% Z-normalization
Z = bsxfun(@minus,corrArr_reg,reg_CN_mean);
Z = bsxfun(@rdivide,Z,reg_CN_std);
% Discretization
discretized_Z = floor(10*Z); % Positive is hyper-connectivity, negative is hypo-connectivity

disp('Z-score done.');

%% Confirm that there is no NaN, Inf, -Inf in your doc
if find(isnan(discretized_Z))
    error('Error: Find NaN in corpus.\n')
end
if find(isinf(discretized_Z))
    error('Error: Find Inf or -Inf in corpus.\n')
end

%% Write into docs for polarLDA
if size(cohort_label,2) == 1
    cohort_label = cohort_label';
end
for dx = cohort_label
    fprintf('---Cohort %d: \n',dx);
    fileID = fopen([output_dir output_name '_dx' num2str(dx) '.dat'], 'w'); % clear contents
    fileID = fopen([output_dir output_name '_dx' num2str(dx) '.dat'], 'a'); % start appending
    for idx1 = 1:no_subjects         
        if dx_info(idx1) == dx
            fprintf('Subject order: %d \n',idx1);
            tc_one_sub = discretized_Z(idx1, :);
            no_terms = sum(tc_one_sub~=0);
            fprintf(fileID, '%i ', no_terms);
            for idx2 = 1:numel(tc_one_sub)
                if tc_one_sub(idx2) ~= 0
                    fprintf(fileID, '%i:%i ', idx2-1, tc_one_sub(idx2));
                end
            end
            fprintf(fileID, '\n');
        end
    end
    fclose(fileID);
end

disp('Writing word counts done.');
