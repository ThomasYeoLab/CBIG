function CBIG_KRDNN_test_krr(setup_file, lambda_set)
% CBIG_KRDNN_test_krr(setup_file, lambda_set)
% 
% This function runs kernel ridge regression algorithm on training and
% testing data.
%
% Inputs:
%   - setup_file
%     Path of the setup file (.txt file).
%     setup_file has the following fields:
%     1. Line 1
%        `data`, the location of your output and intermediate values. You can
%        also change this to any place you want.
%     2. Line 2
%        Path of subject list of UK biobank dataset. It should be a txt file
%        that contains #subject of line (8868 in our paper), while each line
%        is the subject id of 1 subject.
%     3. Line 3
%        Path of CSV file extracted from UK biobank dataset. You can
%        remove the unused measures from the csv to increase the loading
%        speed.
%     4. Line 4
%        Path of the mean rfMRI head motion (data-field 25741) for each
%        subject. It should be a txt file that contains #subject of line,
%        while each line corresponds to the mean rfMRI head motion of a
%        subject. The order should follow the subject order of "subject_list"
%     5. Line 5
%        Path of the functional connectivity matrix. A matrix "corr_mat"
%        is assumed to be saved in this file. "corr_mat" should be a 3-D
%        matrix with dimension of #ROI x #ROI x #subjects. Since it is a
%        connectivity matrix, only the lower-triangular off-diagonal entries
%        will be considered as features because the connectivity matrix is
%        symmetric. The order of "corr_mat" should follow the subject order
%        of "subject_list" 
%     6. Line 6
%        Path of measure list text file
%     7. Line 7
%        Name of the measures lists text file in measure_dir folder. You
%        do not need to change this unless you want to change the measure list
%   - lambda_set
%     1xN vector, where N is number of lambda to search. This is the set of
%     lambda for kernel regression algorithm to search through.
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
assert(isvector(lambda_set), 'lambda_set should be a vector')
[data_dir, ~, ~, ~, ~, ~, measure_set] = CBIG_KRDNN_mics_load_setup(setup_file);
measures_set = {measure_set};
FSM_set = {'corr'};
out_dir = '/test/dis_test';

for i = 1:length(FSM_set)
    curr_FSM_type = FSM_set{i};
    for m = 1:length(lambda_set)
        l = lambda_set(m);
        for j = 1:length(measures_set)
            
            curr_measure = measures_set{j};
            fprintf('FSM: %s  lambda: %d  measure: %s \n', curr_FSM_type,l,curr_measure);
            
            [y_p,y_t,acc_tmp,acc_concat{i,m}] = test_krr_ind(l,curr_measure,curr_FSM_type,out_dir,data_dir);
            if(j == 1)
                acc_all = acc_tmp;
            else
                acc_all = [acc_all acc_tmp];
            end
            acc{i,m} = acc_all;
            avg_acc(i,m) = mean(acc_all);
        end
    end
end
save(fullfile(data_dir, out_dir, 'acc_corr.mat'), 'acc', 'acc_concat', 'avg_acc');
end

function [y_p,y_t,acc,acc_concat] = test_krr_ind(lambda,measure_type,FSM_type,out_dir,data_dir)

% [y_p,y_t,acc,acc_concat] = test_krr_ind(lambda,measure_type,FSM_type,out_dir,data_dir)
% 
% This function performs the kernel regression algorithm with specific
% arguement.
%
% Inputs:
%   - lambda
%     lambda is the regularization parameter for L2 regularization
%   - measure_type
%     type of measures. For UK Biobank dataset, only 1 type and it is used
%     load data.
%   - FSM_type
%     type of kernel. it has three options: 'corr' (Pearson's correlation),
%     'Gaussian' (Gaussian kernel) and 'Exponential' (exponential kernel).
%   - out_dir
%     path of output directory for the inner loop cross validation result
%   - data_dir
%     path of root directory of output and intermediate data
%
% Outputs:
%   - y_p
%     predicted value of behavioral measures y
%   - y_t
%     original value of behavioral measures y
%   - acc
%     accuracy of predicted value of behavioral measures y against original
%     value
%   - acc_concat
%     accuracy of concatenated predicted value of behavioral measures y
%     against original value.

output_dir = fullfile(data_dir, out_dir, 'results', FSM_type);
output_fig_dir = fullfile(data_dir, out_dir, 'figures', FSM_type, ['/lambda' num2str(lambda)]);

%% Extract behaviour measures after regressing out covariates
% after loading, there will be two variable: 
% beta: (1+#covariates) x #measures
% measures_num_regress: #subjects x #measures
load(fullfile(data_dir, 'y', 'dis_y_regress', [measure_type '_y_regress.mat']));

% load train/valid/test split
load(fullfile(data_dir, 'ukbb_subject_split.mat'))
dis_idx = subject_split.index_tra;
cv_idx = subject_split.index_tes;

dis_y_regress = measures_num_regress(dis_idx,:);
cv_y_regress = measures_num_regress(cv_idx,:);
num_subject = size(dis_y_regress,1);
num_subject_cv = size(cv_y_regress,1);

%% Load K 
% K: #subject_dis x #subject_dis
load(fullfile(data_dir, 'FSM_cv', 'dis_test', ['FSM_' FSM_type '.mat']));

% extract behaviour measurements from restriced and unrestriced sheet
measure_csv_file = fullfile(data_dir, 'measures_lists', [measure_type '.txt']);
measure_cell = read_sub_list(measure_csv_file);

if(exist(fullfile(output_dir, [measure_type '_lambda' num2str(lambda) '.mat'])))
    load(fullfile(output_dir, [measure_type '_lambda' num2str(lambda) '.mat']))
    if(exist('acc_fold','var'))
        flag = 0;
    else
        flag = 1;
    end
else
    flag = 1;
end

dis_idx(subject_split.index_val) = [];
cv_idx(subject_split.index_val) = [];

if(flag == 1)
    fold = 1;
    curr_FSM_train = FSM(dis_idx,dis_idx);
    curr_FSM_test = FSM(cv_idx,dis_idx);
    % Training
    % alpha = (K + lambda * I)^-1 * y
    %  Nx1    NxN   1x1    NxN      Nx1
    for i = 1:length(measure_cell)
        curr_y_regress = dis_y_regress(:,i);
        curr_y_real = cv_y_regress(:,i);
        nan_index = isnan(curr_y_regress);
        alpha{fold,i} = (curr_FSM_train(~nan_index,~nan_index) +...
            lambda*eye(sum(~nan_index))) \ curr_y_regress(~nan_index);

        % Test
        y_predict{fold,i} = curr_FSM_test(~isnan(curr_y_real),~nan_index)*alpha{fold,i};

        % Test real
        y_true{fold,i} = curr_y_real(~isnan(curr_y_real));
        acc_fold{fold,i} = corr(y_predict{fold,i},y_true{fold,i});
    end
    mkdir(output_dir);
    save(fullfile(output_dir, [measure_type '_lambda' num2str(lambda) '.mat']),...
        'alpha','y_predict','y_true','acc_fold');
else
    load(fullfile(output_dir, [measure_type '_lambda' num2str(lambda) '.mat']));
end

%% visualization
color_list = jet(10);
measure_cell_name = strrep(measure_cell,'_','\_');
for i = 1:length(measure_cell)
    figure;
    set(gcf, 'Visible', 'off');
    for fold = 1
        hold on;
        scatter(y_predict{fold,i},y_true{fold,i},[],color_list(fold,:));
        
        if (fold == 1)
            y_predict_concat = y_predict{fold,i};
            y_true_concat = y_true{fold,i};
            acc_all = acc_fold{fold,i};
        else
            y_predict_concat = [y_predict_concat; y_predict{fold,i}];
            y_true_concat = [y_true_concat; y_true{fold,i}];
            acc_all = acc_all +  acc_fold{fold,i};
        end
    end
    hold off;
    acc(i) = acc_all;
    acc_concat(i) = corr(y_predict_concat,y_true_concat);
    title(['measure: ' measure_cell_name{i} ' ,r = ' num2str(acc(i)) ' , r\_concat = ' num2str(acc_concat(i))],...
        'FontSize',20,'FontWeight','bold')

    y_p{i} = y_predict_concat;
    y_t{i} = y_true_concat;
    mkdir(output_fig_dir);
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf, fullfile(output_fig_dir, [measure_cell{i} '_lambda' num2str(lambda) '.jpg']));
end
close all;
end

function subj_list = read_sub_list(subject_text_list)
% this function will output a 1xN cell where N is the number of
% subjects in the text_list, each subject will be represented by one
% line in the text file
% NOTE: multiple runs of the same subject will still stay on the same
% line
% Each cell will contain the location of the subject, for e.g.
% '<full_path>/subject1_run1_bold.nii.gz <full_path>/subject1_run2_bold.nii.gz'
    fid = fopen(subject_text_list, 'r');
    i = 0;
    while(1);
        tmp = fgetl(fid);
        if(tmp == -1)
            break
        else
            i = i + 1;
            subj_list{i} = tmp;
        end
    end
    fclose(fid);
end