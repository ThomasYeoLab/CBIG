classdef CBIG_MSHBM_Epilepsy_unit_test < matlab.unittest.TestCase
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)

            % ====================================
            %% Setup
            % ====================================

            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
            proj_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Lim2026_MSHBM_Epilepsy');
            out_dir  = fullfile(proj_dir, 'unit_tests', 'output');
            mkdir(out_dir);

            % ====================================
            %% Run Example Wrapper
            % ====================================

            addpath(fullfile(proj_dir, 'examples'));
            CBIG_MSHBM_Epilepsy_wrapper(out_dir);

            % ====================================
            %% Check Results
            % ====================================

            CBIG_MSHBM_Epilepsy_check_example_results(out_dir);

            % ====================================
            %% Replace Reference (if flag set)
            % ====================================

            ref_dir = fullfile(proj_dir, 'examples', 'results');
            if replace_unittest_flag
                disp('Replacing unit test reference results for CBIG_MSHBM_Epilepsy_unit_test...');
                copyfile(out_dir, ref_dir);
            end

            % ====================================
            %% Cleanup
            % ====================================

            rmdir(out_dir, 's');
            rmpath(fullfile(proj_dir, 'examples'));

        end
    end
end
