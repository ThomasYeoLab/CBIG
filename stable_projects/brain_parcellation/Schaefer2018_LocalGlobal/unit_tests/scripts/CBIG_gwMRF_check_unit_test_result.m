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


ref_dir = '/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/results';


%-------------------------------------------------------------------------%
% check concatenated time matrices
%-------------------------------------------------------------------------%

fprintf('\n[CHECK 1]\t Concatenated time matrices \n');
disp('-------------------------------------------------------------------');

if(~exist([output_dir '/time_data/rh_time_matrix.mat'], 'file'))
	fprintf('[FAILED]\t Time matrix missing! \n');
    
else
   
    % check if lh time matrices are the same
    ref_lh_time = load([ref_dir '/time_data/lh_time_matrix.mat']);
    lh_time = load([output_dir '/time_data/lh_time_matrix.mat']);
    
    if(~isequal(ref_lh_time.lh_time_mat, lh_time.lh_time_mat))
        if(~isequal(size(ref_lh_time.lh_time_mat), size(lh_time.lh_time_mat)))
            fprintf('[FAILED]\t Left hemisphere time matrices are of different size! \n');
        else
            diff_lh_time = max(max(abs(ref_lh_time.lh_time_mat - lh_time.lh_time_mat)));
            fprintf('[FAILED]\t Left hemisphere time matrices are of different value, max_diff = %f \n', diff_lh_time);
        end
    else
        fprintf('[PASSED]\t Left hemisphere time matrices are the same \n');
    end 
    
    % release some memory
    clear ref_lh_time;
    clear lh_time;
     
    % check if rh time matrices are the same
    ref_rh_time = load([ref_dir '/time_data/rh_time_matrix.mat']);
    rh_time = load([output_dir '/time_data/rh_time_matrix.mat']);
     
    if(~isequal(ref_rh_time.rh_time_mat, rh_time.rh_time_mat))
        if(~isequal(size(ref_rh_time.rh_time_mat), size(rh_time.rh_time_mat)))
            fprintf('[FAILED]\t Right hemisphere time matrices are of different size! \n');
        else
            diff_rh_time = max(max(abs(ref_rh_time.rh_time_mat - rh_time.rh_time_mat)));
            fprintf('[FAILED]\t Right hemisphere time matrices are of different value, max_diff = %f \n', diff_rh_time);
        end
    else
         fprintf('[PASSED]\t Right hemisphere time matrices are the same \n');
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

if(~exist([output_dir '/mult_mat/rh_mult_matrix.mat'], 'file'))
	fprintf('[FAILED]\t Premultiplied product matrix missing! \n');
    
else
  
    % check if lh mult matrices are the same
    ref_lh_mult = load([ref_dir '/mult_mat/lh_mult_matrix.mat']);
    lh_mult = load([output_dir '/mult_mat/lh_mult_matrix.mat']);
    
    if(~isequal(ref_lh_mult.cov_mat, lh_mult.cov_mat))
        if(~isequal(size(ref_lh_mult.cov_mat), size(lh_mult.cov_mat)))
            fprintf('[FAILED]\t Left hemisphere mult matrices are of different size! \n');
        else
            diff_lh_mult = max(max(abs(ref_lh_mult.cov_mat-lh_mult.cov_mat))); 
            fprintf('[FAILED]\t Left hemisphere mult matrices are of different value, max_diff = %f \n', diff_lh_mult);
        end
    else
        fprintf('[PASSED]\t Left hemisphere mult matrices are the same \n');
    end 
     
    % release some memory
    clear ref_lh_mult;
    clear lh_mult;
    
    % check if rh mult matrices are the same    
    ref_rh_mult = load([ref_dir '/mult_mat/rh_mult_matrix.mat']);
    rh_mult = load([output_dir '/mult_mat/rh_mult_matrix.mat']);
    
    if(~isequal(ref_rh_mult.cov_mat, rh_mult.cov_mat))
         if(~isequal(size(ref_rh_mult.cov_mat), size(rh_mult.cov_mat)))
            fprintf('[FAILED]\t Right hemisphere mult matrices are of different size! \n');
        else
            diff_rh_mult = max(max(abs(ref_rh_mult.cov_mat-rh_mult.cov_mat)));
            fprintf('[FAILED]\t Right hemisphere mult matrices are of different value, max_diff = %f \n', diff_rh_mult);
        end
    else
         fprintf('[PASSED]\t Right hemisphere mult matrices are the same \n');
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

if(~exist([output_dir '/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_7_seed_2.mat'], 'file'))
	fprintf('[FAILED]\t Parcellation result missing! \n'); % there should be two seeds in 'output_dir/clustering'
    
else
    % start comparing seed 1
    ref_seed_1 = load([ref_dir '/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_7_seed_1.mat']);
    seed_1 = load([output_dir '/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_7_seed_1.mat']);
   
    [output, assign, cost_lh, dice_overlap] = CBIG_HungarianClusterMatch(ref_seed_1.results.lh_label, seed_1.results.lh_label', 0);
    [output, assign, cost_rh, dice_overlap] = CBIG_HungarianClusterMatch(ref_seed_1.results.rh_label, seed_1.results.rh_label', 0);
    
    if(abs(cost_lh) ~= 37476)
        fprintf('[FAILED]\t Left hemisphere parcellations of seed 1 are diffrent, overlap_percentage = %f \n', abs(cost_lh)/37476);
    else
        fprintf('[PASSED]\t Left hemisphere parcellations of seed 1 are the same \n');
    end
    
    if(abs(cost_rh) ~= 37471)
        fprintf('[FAILED]\t Right hemisphere parcellations of seed 1 are diffrent, overlap_percentage = %f \n', abs(cost_rh)/37471);
    else
        fprintf('[PASSED]\t Right hemisphere parcellations of seed 1 are the same \n');
    end
    
    clear cost_lh;
    clear cost_rh;

    % start comparing seed 2
    ref_seed_2 = load([ref_dir '/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_7_seed_2.mat']);
    seed_2 = load([output_dir '/clustering/Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_7_seed_2.mat']);

    [output, assign, cost_lh, dice_overlap] = CBIG_HungarianClusterMatch(ref_seed_2.results.lh_label, seed_2.results.lh_label', 0);
    [output, assign, cost_rh, dice_overlap] = CBIG_HungarianClusterMatch(ref_seed_2.results.rh_label, seed_2.results.rh_label', 0);
    
    if(abs(cost_lh) ~= 37476)
        fprintf('[FAILED]\t Left hemisphere parcellations of seed 2 are diffrent, overlap_percentage = %f \n', abs(cost_lh)/37476);
    else
        fprintf('[PASSED]\t Left hemisphere parcellations of seed 2 are the same \n');
    end
    
    if(abs(cost_rh) ~= 37471)
        fprintf('[FAILED]\t Right hemisphere parcellations of seed 2 are diffrent, overlap_percentage = %f \n', abs(cost_rh)/37471);
    else
        fprintf('[PASSED]\t Right hemisphere parcellations of seed 2 are the same \n');
    end
    
end

end
