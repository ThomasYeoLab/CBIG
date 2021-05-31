classdef CBIG_ArealMSHBM_unit_test < matlab.unittest.TestCase
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2022_ArealMSHBM', 'unit_tests');
            out_dir = fullfile(cur_dir, 'output');
            mkdir(out_dir)
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2022_ArealMSHBM', 'examples'));
            
            CBIG_ArealMSHBM_example_wrapper(out_dir);
            
            % check example results
            CBIG_ArealMSHBM_check_example_results(out_dir, 1);
            
            % remove the output directory
            rmdir(out_dir, 's')

            % remove path
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2022_ArealMSHBM', 'examples'));
        end
    end
end