function CBIG_KRDNN_get_ukbb_krr_result(setup_file, lambda_set)
% CBIG_KRDNN_get_ukbb_krr_result(setup_file, lambda_set)
% 
% This function get the result of kernel ridge regression based on the result
% of validation result and testing result.
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
    type = 'corr';
    FSM_set = {type};

    innerloop_acc = load(fullfile(data_dir, 'valid', 'dis_valid', ['acc_' type '.mat']));
    test_acc = load(fullfile(data_dir, 'test', 'dis_test', ['acc_' type '.mat']));
    innerloop_acc_out = rearrange_acc(innerloop_acc.acc);
    test_acc_out = rearrange_acc(test_acc.acc);
    for i = 1:4  % 4 is number of behavioral measures 
        curr_innerloop_acc = innerloop_acc_out(:,:,i);
        curr_test_acc = test_acc_out(:,:,i);
        [~,I] = max(curr_innerloop_acc(:));
        [FSM_idx,lambda_idx] = ind2sub(size(curr_innerloop_acc),I);
        opt_FSM = FSM_set{FSM_idx};
        opt_lambda = lambda_set(lambda_idx);
        
        load(fullfile(data_dir, 'test', 'dis_test', 'results', opt_FSM,...
            [measure_set '_lambda' num2str(opt_lambda) '.mat']));
        curr_y_predict{i} = y_predict{i};
        curr_y_true{i} = y_true{i};

        if i == 1
            sex_accuracy = accuracy_for_sex(data_dir,...
                fullfile('results', opt_FSM, [measure_set '_lambda' num2str(opt_lambda) '.mat']));
        end
        if i == 2
            age_mae = mean(abs(y_predict{i} - y_true{i}));
        end
        optimal_acc(i) = curr_test_acc(I);
    end

    optimal_corr = optimal_acc(2:end);
    save(fullfile(data_dir, 'final_result.mat'), 'optimal_corr', ...
        'sex_accuracy', 'age_mae', 'curr_y_predict', 'curr_y_true')
end

function result = accuracy_for_sex(data_dir, file_name)
% result = accuracy_for_sex(data_dir, file_name)
% 
% This function generate the accuracy value for sex prediction
%
% Inputs:
%   - data_dir
%     path of root directory of output and intermediate data
%   - file_name
%     name of mat file which saved the original and predicted sex value
%
% Outputs:
%   - result
%     accuracy of sex prediction for the data in file 'file_name'

    load(fullfile(data_dir ,'valid', 'dis_valid', file_name))
    y_true = y_true{1};
    y_predict = y_predict{1};
    male_female = unique(y_true);

    assert(size(male_female,1) == 2, 'Validation set: y_ture does not have 2 unique value')

    low = min(male_female);
    high = max(male_female);

    threshold = 0;
    record = 0;
    for i = low:0.001:high
        temp = y_predict;
        temp(temp>i) = high;
        temp(temp<=i) = low;
        score = 1 - floor(sum(abs(y_true-temp))) / size(y_true, 1);
        if score >= record
            record = score;
            threshold = i;
        end
    end

    load(fullfile(data_dir ,'test', 'dis_test', file_name))

    y_true = y_true{1};
    y_predict = y_predict{1};
    male_female = unique(y_true);

    assert(size(male_female,1) == 2, 'Test set: y_ture does not have 2 unique value')
    low = min(male_female);
    high = max(male_female);
    temp = y_predict;
    temp(temp>threshold) = high;
    temp(temp<=threshold) = low;
    result = 1 - floor(sum(abs(y_true-temp))) / size(y_true, 1);
end

function acc_out = rearrange_acc(acc)
% acc_out = rearrange_acc(acc)
% 
% This function rearrange input acc (cell) to matrix
%
% Inputs:
%   - acc
%     input acc in cell format
%
% Outputs:
%   - acc_out
%     matrix version of input acc

    for i = 1:size(acc,1)
        for j = 1:size(acc,2)
            acc_each = acc{i,j};
            for m = 1:length(acc_each)
                acc_out(i,j,m) = acc_each(m);
            end
        end
    end
end