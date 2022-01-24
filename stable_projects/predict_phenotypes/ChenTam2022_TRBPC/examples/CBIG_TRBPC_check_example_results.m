function CBIG_TRBPC_check_example_results(results_dir)

% function CBIG_TRBPC_check_example_results(results_dir, ref_dir)
%
% This function checks if example output are the same as the reference
% output and throw an error when results are difference from reference
%
% Inputs:
%   -results_dir
%   example results directory from CBIG_TRBPC_example_wrapper.m
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

example_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes','ChenTam2022_TRBPC','examples');
ref_dir = fullfile(example_dir, 'ref_output');

% check single-kernel KRR
ref_file = fullfile(ref_dir,'singleKRR','final_result_all_score.mat');
result_file = fullfile(results_dir,'singleKRR','final_result_all_score.mat');
check_KRR_results(ref_file, result_file);

% check multi-kernel KRR
ref_file = fullfile(ref_dir,'multiKRR','final_result_all_score.mat');
result_file = fullfile(results_dir,'multiKRR','final_result_all_score.mat');
check_KRR_results(ref_file, result_file);

% check LRR
ref_file = fullfile(ref_dir,'LRR','final_result_all_score.mat');
result_file = fullfile(results_dir,'LRR','final_result_all_score.mat');
check_LRR_results(ref_file, result_file);

% check single-kernel KRR PFM
check_PFM_results(fullfile(ref_dir,'PFM','singleKRR'), fullfile(results_dir, 'PFM', 'singleKRR'));

% check multi-kernel KRR PFM
check_PFM_results(fullfile(ref_dir,'PFM','multiKRR'), fullfile(results_dir, 'PFM', 'multiKRR'));

% check LRR PFM
check_PFM_results(fullfile(ref_dir,'PFM','LRR'), fullfile(results_dir, 'PFM', 'LRR'));

end

function check_KRR_results(result_file, ref_file)

ref = load(ref_file);
test = load(result_file);
fields = fieldnames(ref);

for i = 1:length(fields)
    if  isequal(fields{i}, 'optimal_y_pred') || isequal(fields{i}, 'optimal_hyp')
        curr_test = test.(fields{i});
        curr_ref = ref.(fields{i});
        for j = 1:length(curr_ref)
            test_fold = curr_test{j};
            ref_fold = curr_ref{j};
            for k = 1:length(ref_fold)
                assert(isequal(size(test_fold{k}),size(ref_fold{k})), ...
                    [fields{i} ': size is different from reference']);
                assert(max(abs((test_fold{k} - ref_fold{k}))) < 1e-6, ...
                    [fields{i} ': different from reference']);
            end
        end
        
    elseif isequal(fields{i}, 'optimal_stats')
        stats_ref = ref.optimal_stats;
        stats_test = test.optimal_stats;
        stats_names = fieldnames(stats_ref);
        for j = 1:length(stats_names)
            curr_stats_ref = stats_ref.(stats_names{j});
            curr_stats_test = stats_test.(stats_names{j});
            assert(max(abs(curr_stats_ref(:) - curr_stats_test(:))) < 1e-6, ...
                'optimal stats are different from reference result.');
        end
        
    elseif isequal(fields{i}, 'optimal_kernel')
        assert(isequaln(test.optimal_kernel, ref.optimal_kernel),'optimal kernel is different')
        
    elseif isequal(fields{i}, 'optimal_threshold')
        if ~all(isnan(ref.optimal_threshold(:)))
            mask = ~isnan(ref.optimal_threshold);
            thr_ref = ref.optimal_threshold(mask);
            thr_test = test.optimal_threshold(mask);
            assert(max(abs((thr_ref(:) - thr_test(:)))) < 1e-6, 'optimal thresholds different');
        end
    elseif isequal(fields{i},'y_pred_train')
        curr_ref = ref.(fields{i});
        curr_test = test.(fields{i});
        for j = 1:length(curr_ref)
            assert(max(abs((curr_test{j}(:) - curr_ref{j}(:)))) < 1e-6, ...
            sprintf('field %s is different from reference result.', fields{i}));
        end
    else
        curr_ref = ref.(fields{i});
        curr_test = test.(fields{i});

        assert(isequal(size(curr_test),size(curr_ref)), ...
            sprintf('field %s is of wrong size.', fields{i}));
        assert(max(abs((curr_test(:) - curr_ref(:)))) < 1e-6, ...
            sprintf('field %s is different from reference result.', fields{i}));
    end
end

disp('All results replicated')

end

function check_LRR_results(result_file, ref_file)

ref = load(ref_file);
test = load(result_file);
fields = fieldnames(ref);

for i = 1:length(fields)
    if  isequal(fields{i}, 'optimal_y_pred')
        for j = 1:size(ref.optimal_y_pred,1)
            for k = 1:size(ref.optimal_y_pred,2)
                curr_ref = ref.optimal_y_pred{j,k};
                curr_test = test.optimal_y_pred{j,k};
                assert(max(abs((curr_test(:) - curr_ref(:)))) < 1e-6, ...
                    [fields{i} ': different from referecne']);
            end
        end
    elseif isequal(fields{i}, 'optimal_stats')
        stats_ref = ref.optimal_stats;
        stats_test = test.optimal_stats;
        stats_names = fieldnames(stats_ref);
        for j = 1:length(stats_names)
            curr_stats_ref = stats_ref.(stats_names{j});
            curr_stats_test = stats_test.(stats_names{j});
            assert(max(abs(curr_stats_ref(:) - curr_stats_test(:))) < 1e-6, ...
                'optimal stasts are different from reference result.');
        end
    else
        curr_ref = ref.(fields{i});
        curr_test = test.(fields{i});

        assert(isequal(size(curr_test),size(curr_ref)), ...
            sprintf('field %s is of wrong size.', fields{i}));
        assert(max(abs((curr_test(:) - curr_ref(:)))) < 1e-6, ...
            sprintf('field %s is different from reference result.', fields{i}));
    end
end
end

function check_PFM_results(results_dir, ref_dir)

for i = 1:2
    test = load(fullfile(results_dir, ['PFM_score' num2str(i) '_all_folds.mat']));
    ref = load(fullfile(ref_dir, ['PFM_score' num2str(i) '_all_folds.mat']));
    assert(max(abs((test.PFM_all_folds(:) - ref.PFM_all_folds(:)))) < 1e-6, ...
        'PFM different');
end

end