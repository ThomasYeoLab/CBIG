classdef CBIG_EIAD_unit_test < matlab.unittest.TestCase
    % Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag.txt'));

            proj_dir   = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'fMRI_dynamics', 'Zhang2026_EIAD');
            script_dir = fullfile(proj_dir, 'examples', 'script');
            out_dir    = fullfile(proj_dir, 'examples', 'output');
            ref_dir    = fullfile(proj_dir, 'examples', 'reference_output');
            out_file   = fullfile(out_dir, 'example_output.mat');
            ref_file   = fullfile(ref_dir, 'reference_output.mat');

            cur_dir = pwd;
            addpath(script_dir);
            cd(script_dir);

            % Run the example WITHOUT its own comparison (run_check = false);
            % this only produces examples/output/example_output.mat.
            CBIG_EIAD_example(false);

            cd(cur_dir);

            % load the results of this run
            example_output = load(out_file);
            output_struct  = example_output.output_struct;

            if(replace_unittest_flag)
                % replace the reference output with the results of this run
                disp('Replacing unit test reference results for CBIG_EIAD_unit_test...');
                copyfile(out_file, ref_file);
            else
                % compare against the reference output with the dedicated checker,
                % which throws an error if any field differs beyond the tolerance
                try
                    CBIG_EIAD_check_example_results(output_struct, ref_file);
                    check_error = '';
                catch ME
                    check_error = ME.message;
                end
                testCase.verifyEmpty(check_error, check_error);
            end

            % remove the output directory
            if(exist(out_dir, 'dir'))
                rmdir(out_dir, 's')
            end

            % remove path
            rmpath(script_dir);
        end
    end
end
