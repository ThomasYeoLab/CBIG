classdef CBIG_GradPar_unit_test < matlab.unittest.TestCase
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));

            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'Kong2023_GradPar', 'unit_tests');
            out_dir = fullfile(cur_dir, 'output');
            mkdir(out_dir)
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'Kong2023_GradPar', 'examples'));
            
            CBIG_GradPar_example_LRR_frac_wrapper(out_dir);
            CBIG_GradPar_example_KRR_wrapper(out_dir);
            
            % check example results
            CBIG_GradPar_check_example_results(out_dir);
            
            % replace reference output if flag is 1
            ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'Kong2023_GradPar', 'examples', 'ref_results');
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_GradPar_unit_test...');
                copyfile(out_dir, ref_dir);
            end

            % remove the output directory
            rmdir(out_dir, 's')

            % remove path
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'Kong2023_GradPar', 'examples'));
        end
    end
end