classdef CBIG_TRBPC_unit_test < matlab.unittest.TestCase
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            project_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'ChenTam2022_TRBPC');
            example_dir = fullfile(project_dir, 'examples');
            outdir = fullfile(project_dir, 'unit_tests', 'output');
            if exist(outdir,'dir')
                rmdir(outdir, 's')
            end
            mkdir(outdir)
            
            % run the example
            addpath(example_dir)
            CBIG_TRBPC_example_wrapper(outdir)
            
            % compare the results
            CBIG_TRBPC_check_example_results(outdir)

            % remove the output directory
            rmdir(outdir, 's')
            rmpath(example_dir)
        end
    end
end