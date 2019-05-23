function [] = CBIG_gwMRF_check_unit_test_result(output_dir)

% [] = CBIG_gwMRF_check_unit_test_result(output_dir)
%
% This function compares the output of unit test with the ref output
% on our HPC server to make sure the codes work fine.
%
% '[FAILED]' message indicates something wrong with the codes, the
% unit test is successful only if all messages are '[PASSED]'.
%   
% Input:
%   - output_dir = the folder containing the output of your unit test
%
% Example:
%   CBIG_gwMRF_check_unit_test_result('~/storage/Temporary/CBIG_gwMRF_unit_test');
%
% Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


unit_test_data_dir = fullfile('/mnt', 'eql', 'yeo1', 'CBIG_private_data', 'unit_tests');
ref_dir = fullfile(unit_test_data_dir, 'stable_projects', ...
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
file_seed1 = [prefix, '_seed_1_lh_reduction_iteration_1.mat'];
file_seed2 = [prefix, '_seed_2_rh_reduction_iteration_1.mat'];

if(~exist(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed1), 'file') ...
        || ~exist(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed2), 'file'))
	fprintf('[FAILED]\t Parcellation result missing! \n');
else
    % start comparing seed 1
    ref_seed_1 = load(fullfile(ref_dir, 'clustering', 'inbetween_results', file_seed1));
    seed_1 = load(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed1));
   
    [~, ~, cost, ~] = CBIG_HungarianClusterMatch(ref_seed_1.current_label,...
        seed_1.current_label', 0);
    
    if(abs(cost) ~= 37476)
        fprintf('[FAILED]\t LH intermediate parcellations of seed 1 are diffrent, overlap_percentage = %f \n', ...
            abs(cost_lh)/37476);
    else
        fprintf('[PASSED]\t LH intermediate parcellations of seed 1 are the same \n');
    end
    
    clear cost;

    % start comparing seed 2
    ref_seed_2 = load(fullfile(ref_dir, 'clustering', 'inbetween_results', file_seed2));
    seed_2 = load(fullfile(output_dir, 'clustering', 'inbetween_results', file_seed2));

    [~, ~, cost, ~] = CBIG_HungarianClusterMatch(ref_seed_2.current_label, ...
        seed_2.current_label', 0);
    
    if(abs(cost) ~= 37471)
        fprintf('[FAILED]\t RH intermediate parcellations of seed 2 are diffrent, overlap_percentage = %f \n', ...
            abs(cost_rh)/37471);
    else
        fprintf('[PASSED]\t RH intermediate parcellations of seed 2 are the same \n');
    end
    
end

end
