classdef CBIG_MM_unit_test < matlab.unittest.TestCase
% Written by He Tong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_MM(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2022_MM', 'unit_tests');

            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2022_MM', 'examples'));
            
            CBIG_MM_example();

            % get output folder
            out_dir = fullfile(cur_dir, 'output_KRR');
            
            CBIG_MM_check_example_results(out_dir, false);
            
            % remove the output directory and path
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2022_MM', 'examples'));
        end
    end
end