classdef CBIG_IndCBM_unit_test < matlab.unittest.TestCase
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Xue2021_IndCerebellum', 'unit_tests');
            output_dir = fullfile(cur_dir, 'output');
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Xue2021_IndCerebellum', 'examples'));
            
            CBIG_IndCBM_example_wrapper(output_dir);
            
            % check example results
            CBIG_IndCBM_check_example_results(output_dir);
            
            % remove the output directory
            rmdir(output_dir, 's')

            % remove path
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Xue2021_IndCerebellum', 'examples'));
        end
    end
end
