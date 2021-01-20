classdef CBIG_KRDNN_unit_test < matlab.unittest.TestCase
% Written by He Tong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_KR_TVT(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'unit_tests');

            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
            
            CBIG_KRDNN_KR_TVT_example();

            % get output folder
            out_dir = fullfile(cur_dir, 'output', 'output_kr_tvt_example');
            
            CBIG_KRDNN_check_example_results(out_dir, false);
            
            % remove the output directory and path
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
        end

        function test_KR_CV(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'unit_tests');

            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
            
            CBIG_KRDNN_KR_CV_example();

            % get output folder
            out_dir = fullfile(cur_dir, 'output', 'output_kr_cv_example');
            
            CBIG_KRDNN_check_example_results(out_dir, true);
            
            % remove the output directory and path
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
            rmdir(fullfile(cur_dir, 'output'))
        end
    end
end