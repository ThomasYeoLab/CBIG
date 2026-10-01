classdef CBIG_LBC_unit_test < matlab.unittest.TestCase
% Run unit test for Xie2025_LBC project. Tests whether KRR prediction
% accuracies is consistent with results.
% Written by Yapei Xie and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % run through KRR prediction results using simulated data
            % and check whether results replicate

            % set directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            project_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Xie2025_LBC');
        
            example_dir = fullfile(project_dir, 'examples');
            outdir = fullfile(project_dir, 'unit_tests', 'output');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag.txt'));
            
            addpath(example_dir)

            % create output folder
            if exist(outdir,'dir')
                rmdir(outdir, 's')
            end
            mkdir(outdir)
            
            % run the example
            CBIG_LBC_KRR_example_wrapper(outdir)
            
            % replace unit test results
            if replace_unit_test
                fprintf('Replacing unit test results... \n')
                copyfile(fullfile(outdir), fullfile(example_dir, 'ref_results'))
            end

            % compare the results
            CBIG_LBC_KRR_check_example_results(outdir)

            % remove the output directory
            rmdir(outdir, 's')
            rmpath(example_dir)
        end
    end
end
