function test_array = CBIG_MMP_check_example_results(results_dir, level)

% function CBIG_MMP_check_example_results(results_dir)
%
% This function checks if example output are the same as the reference
% output and throw an error when results are difference from reference
%
% Inputs:
%   -results_dir
%   example results directory from CBIG_TRBPC_example_wrapper.m
%
%   -level
%   Either 1 or 2. Indicates whether to check first-level or second-level results
%
% Outputs:
%   -test_array
%   Array indicating which test did not match the reference result. The array is a 
%   vector of length 6 when checking first level results and vector of length 2 when
%   checking second level results.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

example_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes','Ooi2022_MMP','examples');
ref_dir = fullfile(example_dir, 'ref_output');

% check first level results
if level == 1
    fprintf('Checking first-level results... \n')
    test_array = zeros(6,1);
    reg_types = {'KRR' 'LRR' 'Elasticnet'};
    % read output results
    for i = 1:length(reg_types)
        output = load(fullfile(results_dir, strcat(reg_types{i},'_corr_results.mat')));
        ref = load(fullfile(ref_dir, strcat(reg_types{i},'_corr_results.mat')));
        % check t1 results
        t1_diff = abs(sum(ref.results.t1 - output.results.t1));
        % check fmri results
        fmri_diff = abs(sum(ref.results.fmri - output.results.fmri));
        
        % print difference if more than threshold
        if t1_diff > 1e-8
            fprintf('Sum of %s t1 absolute difference between output and reference is: %f \n', ...
                reg_types{i}, t1_diff)
        else
            test_array((i-1)*2+1,1) = 1;
        end
        if fmri_diff > 1e-8
            fprintf('Sum of %s fmri absolute difference between output and reference is: %f \n', ...
                reg_types{i}, fmri_diff)
        else
            test_array((i-1)*2+2,1) = 1;
        end
    end

% check second level results
elseif level == 2
    fprintf('Checking second-level results... \n')
    test_array = zeros(2,1);
    reg_types = {'multiKRR' 'stacking'};
    output = load(fullfile(results_dir, strcat('combined_models_corr_results.mat')));
    ref = load(fullfile(ref_dir, strcat('combined_models_corr_results.mat')));
    % read output results
    for i = 1:length(reg_types)
        % check combined model results
        model_diff = abs(sum(ref.results.combined(:,:,i+1) - output.results.combined(:,:,i+1)));
        
        % print difference if more than threshold
        if model_diff > 1e-8
            fprintf('Sum of %s (combining cv and rs) absolute difference between output and reference is: %f \n', ...
                reg_types{i}, model_diff)
        else
            test_array(i,1) = 1;
        end
    end
end