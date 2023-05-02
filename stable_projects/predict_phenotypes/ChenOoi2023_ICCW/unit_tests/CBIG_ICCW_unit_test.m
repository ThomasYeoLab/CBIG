classdef CBIG_ICCW_unit_test < matlab.unittest.TestCase
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % run through example and check whether results replicate
            % set directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            project_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'ChenOoi2023_ICCW');
            example_dir = fullfile(project_dir, 'examples');
            outdir = fullfile(project_dir, 'unit_tests', 'output');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            addpath(genpath(example_dir))

            % create output folder
            if exist(outdir,'dir')
                rmdir(outdir, 's')
            end
            mkdir(outdir)
            
            % run the example
            CBIG_ICCW_KRR_example_wrapper(outdir)
            
            % replace unit test results
            if replace_unit_test
                fprintf('Replacing unit test results... \n')
                copyfile(fullfile(outdir), fullfile(example_dir, 'ref_results'))
            end

            % compare the results
            CBIG_ICCW_check_example_results(outdir)

            % remove the output directory
            rmdir(outdir, 's')
            rmpath(genpath(example_dir))
        end
    end
end
