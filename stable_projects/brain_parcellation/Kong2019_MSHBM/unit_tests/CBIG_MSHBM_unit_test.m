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
            
            % check example results
            CBIG_MSHBM_check_example_results(out_dir);
            
            % remove the output directory
            rmdir(out_dir, 's')

            % remove path
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2019_MSHBM', 'examples'));
        end
    end
end