classdef CBIG_MSHBM_unit_test < matlab.unittest.TestCase
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2019_MSHBM', 'unit_tests');
            out_dir = fullfile(cur_dir, 'output');
            mkdir(out_dir)
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2019_MSHBM', 'examples'));
            
            [~, lh_labels, rh_labels] = CBIG_MSHBM_example_wrapper(out_dir);
            
            % check the results
            ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2019_MSHBM', 'examples', 'results');
            
            % compare estimated group priors
            ref_Params = load(fullfile(ref_dir, 'estimate_group_priors', 'priors', 'Params_Final.mat'));
            test_Params = load(fullfile(out_dir, 'estimate_group_priors', 'priors', 'Params_Final.mat'));
          
            diff_mu = max(max(abs(ref_Params.Params.mu - test_Params.Params.mu)));
            assert(diff_mu < 1e-6, sprintf('maximum Params.mu difference: %f',diff_mu));
            
            % compare estimated individual parcellation
            ref_labels = load(fullfile(ref_dir, 'generate_individual_parcellations',...
                'ind_parcellation', 'Ind_parcellation_MSHBM_sub1_w100_MRF50.mat'));
            
            diff_labels = max(abs([lh_labels; rh_labels] - [ref_labels.lh_labels; ref_labels.rh_labels]));
            assert(diff_labels == 0, sprintf('maximum labels difference: %f',diff_labels));
            
            % remove the output directory
            rmdir(out_dir, 's')
        end
    end
end