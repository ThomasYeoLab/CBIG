function [] = CBIG_gwMRF_check_unit_test_result(output_dir)

% [] = CBIG_gwMRF_check_unit_test_result(output_dir)
%
% This function compares the output of unit test with the ref output
% on our HPC server to make sure the codes work fine.
%
% '[FAILED]' message indicates something wrong with the codes, the
% unit test is considered successful only if there is no [FAILED] and the test is [DONE] in the log file.
%   
% Input:
%   - output_dir = the folder containing the output of your unit test
%
% Example:
%   CBIG_gwMRF_check_unit_test_result('~/storage/Temporary/CBIG_gwMRF_unit_test');
%
% Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% check if replacing example results
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));

ref_dir = fullfile(getenv('CBIG_TESTDATA_DIR'), 'stable_projects', ...
    'brain_parcellation', 'Schaefer2018_LocalGlobal', 'results');

%-------------------------------------------------------------------------%
% check concatenated time matrices
%-------------------------------------------------------------------------%

fprintf('\n[CHECK 1]\t Concatenated time matrices \n');
disp('-------------------------------------------------------------------');

if(~exist(fullfile(output_dir, 'time_data', 'rh_time_matrix.mat'), 'file'))
    fprintf('[FAILED]\t Time matrix missing! \n');
else
    % check if lh time matrices are the same
    ref_lh_time = load(fullfile(ref_dir, 'time_data', 'lh_time_matrix.mat'));
    lh_time = load(fullfile(output_dir, 'time_data', 'lh_time_matrix.mat'));
    
    if(~isequal(ref_lh_time.lh_time_mat, lh_time.lh_time_mat))
        if(~isequal(size(ref_lh_time.lh_time_mat), size(lh_time.lh_time_mat)))
            fprintf('[FAILED]\t LH time matrices are of different size! \n');
        else
            diff_lh_time = max(max(abs(ref_lh_time.lh_time_mat - lh_time.lh_time_mat)));
            fprintf('[FAILED]\t LH time matrices are different, max_diff = %f \n',...
                diff_lh_time);
        end
    else
        fprintf('[PASSED]\t LH time matrices are the same \n');
    end 
    
    if(replace_unittest_flag)
        disp('Replacing reference results for lh concatenated time matrices...');
        copyfile(fullfile(output_dir, 'time_data', 'lh_time_matrix.mat'),...
         fullfile(ref_dir, 'time_data', 'lh_time_matrix.mat'));
    end
    % release some memory
    clear ref_lh_time;
    clear lh_time;
     
    % check if rh time matrices are the same
    ref_rh_time = load(fullfile(ref_dir, 'time_data', 'rh_time_matrix.mat'));
    rh_time = load(fullfile(output_dir, 'time_data', 'rh_time_matrix.mat'));
     
    if(~isequal(ref_rh_time.rh_time_mat, rh_time.rh_time_mat))
        if(~isequal(size(ref_rh_time.rh_time_mat), size(rh_time.rh_time_mat)))
            fprintf('[FAILED]\t RH time matrices are of different size! \n');
        else
            diff_rh_time = max(max(abs(ref_rh_time.rh_time_mat - rh_time.rh_time_mat)));
            fprintf('[FAILED]\t RH time matrices are different, max_diff = %f \n', ...
                diff_rh_time);
        end
    else
         fprintf('[PASSED]\t RH time matrices are the same \n');
    end

        
    if(replace_unittest_flag)
        disp('Replacing reference results for rh concatenated time matrices...');
        copyfile(fullfile(output_dir, 'time_data', 'rh_time_matrix.mat'),...
         fullfile(ref_dir, 'time_data', 'rh_time_matrix.mat'));
    end
    % release some memory
    clear ref_rh_time;
    clear rh_time;
    
end

%-------------------------------------------------------------------------%
% check premultiplied product matrices
%-------------------------------------------------------------------------%

fprintf('\n[CHECK 2]\t Premultiplied product matrices \n');
disp('-------------------------------------------------------------------');

if(~exist(fullfile(output_dir, 'mult_mat', 'rh_mult_matrix.mat'), 'file'))
    fprintf('[FAILED]\t Premultiplied product matrix missing! \n');
    
else
    % check if lh mult matrices are the same
    ref_lh_mult = load(fullfile(ref_dir, 'mult_mat', 'lh_mult_matrix.mat'));
    lh_mult = load(fullfile(output_dir, 'mult_mat', 'lh_mult_matrix.mat'));
    
    if(~isequal(ref_lh_mult.cov_mat, lh_mult.cov_mat))
        if(~isequal(size(ref_lh_mult.cov_mat), size(lh_mult.cov_mat)))
            fprintf('[FAILED]\t LH mult matrices are of different size! \n');
        elseif(max(max(abs(ref_lh_mult.cov_mat-lh_mult.cov_mat))) > 1e-6)
            diff_lh_mult = max(max(abs(ref_lh_mult.cov_mat-lh_mult.cov_mat)));
            fprintf('[FAILED]\t LH mult matrices are different, max_diff = %f \n', ...
                diff_lh_mult);
        else
            fprintf('[PASSED]\t LH mult matrices are the same \n');
        end
    else
        fprintf('[PASSED]\t LH mult matrices are the same \n');
    end 
     
    if(replace_unittest_flag)
        disp('Replacing reference results for lh premultiplied product matrices...');
        copyfile(load(fullfile(output_dir, 'mult_mat', 'lh_mult_matrix.mat')),...
         fullfile(ref_dir, 'mult_mat', 'lh_mult_matrix.mat'));
    end
    % release some memory
    clear ref_lh_mult;
    clear lh_mult;
    
    % check if rh mult matrices are the same    
    ref_rh_mult = load(fullfile(ref_dir, 'mult_mat', 'rh_mult_matrix.mat'));
    rh_mult = load(fullfile(output_dir, 'mult_mat', 'rh_mult_matrix.mat'));
    
    if(~isequal(ref_rh_mult.cov_mat, rh_mult.cov_mat))
        if(~isequal(size(ref_rh_mult.cov_mat), size(rh_mult.cov_mat)))
            fprintf('[FAILED]\t RH mult matrices are of different size! \n');
        elseif(max(max(abs(ref_rh_mult.cov_mat-rh_mult.cov_mat))) > 1e-6)
            diff_rh_mult = max(max(abs(ref_rh_mult.cov_mat-rh_mult.cov_mat)));
            fprintf('[FAILED]\t RH mult matrices are different, max_diff = %f \n', ...
                diff_rh_mult);
        else
            fprintf('[PASSED]\t RH mult matrices are the same \n');
        end
    else
         fprintf('[PASSED]\t RH mult matrices are the same \n');
    end
    
    if(replace_unittest_flag)
        disp('Replacing reference results for rh premultiplied product matrices...');
        copyfile(load(fullfile(output_dir, 'mult_mat', 'rh_mult_matrix.mat')),...
         fullfile(ref_dir, 'mult_mat', 'rh_mult_matrix.mat'));
    end
    % release some memory
    clear ref_rh_mult;
    clear rh_mult;
end

%-------------------------------------------------------------------------%
% check parcellation results
%-------------------------------------------------------------------------%

fprintf('\n[CHECK 3]\t Parcellation results \n');
disp('-------------------------------------------------------------------');

prefix = 'Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_2';
file_seed2 = [prefix, '_seed_2_rh_reduction_iteration_1.mat'];

if(~exist(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed2), 'file'))
    fprintf('[FAILED]\t Parcellation result missing! \n');
else
    % start comparing seed 2
    ref_seed_2 = load(fullfile(ref_dir, 'clustering', 'inbetween_results', file_seed2));
    seed_2 = load(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed2));

    [~, ~, cost, ~] = CBIG_HungarianClusterMatch(ref_seed_2.current_label, ...
        seed_2.current_label', 0);
    
    labels_diff = 37471 - abs(cost);
    if(labels_diff == 0)
        fprintf('[PASSED]\t RH intermediate parcellations of seed 2 are the same \n');
    elseif(labels_diff/37471 < 0.0005) % Less than 0.05% of total voxels are different.
        fprintf(['[PASSED]\t' num2str(labels_diff) ' labels are different.\n']);
        warning('The small difference might be caused by different running environments.');
    else % Too many voxels are different.
        fprintf('[FAILED]\t RH intermediate parcellations of seed 2 are different, overlap_percentage = %f \n', ...
            abs(cost_rh)/37471);
    end
end

if(replace_unittest_flag)
    disp('Replacing reference results for parcellation labels...');
    copyfile(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed2),...
    fullfile(ref_dir, 'clustering', 'inbetween_results', file_seed2));
end

fprintf('[DONE]\t Unit test is done. \n');

end
