classdef CBIG_MMP_unit_test < matlab.unittest.TestCase
% Placeholder for unit test
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example_first_level(testCase)
            % check whether first level models match: KRR, LRR and
            % Elasticnet
            
            % set directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            project_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Ooi2022_MMP');
            example_dir = fullfile(project_dir, 'examples');
            outdir = fullfile(project_dir, 'unit_tests', 'output');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(outdir,'dir')
                rmdir(outdir, 's')
            end
            mkdir(outdir)
            
            % run the example
            command = [example_dir, ...
                '/CBIG_MMP_HCP_example_singlefeature_regression_wrapper.sh ' ...
                outdir];
            [status, cmdout] = system(command);

            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N Ooi2022_Example | wc -l'";
            [status, cmdout] = system(check_command);
            while str2num(cmdout) >= 1
                fprintf("Waiting for job(s) to finish... \n")
                pause(30)
                [status, cmdout] = system(check_command);
            end
            
            % compare the results
            addpath(genpath(example_dir))
            reg_types = {'KRR' 'LRR' 'Elasticnet'};
            for i = 1:length(reg_types)
                CBIG_MMP_HCP_collate_result_example_wrapper( ...
                    fullfile(outdir, 'output'), fullfile(outdir, 'results'), ...
                    reg_types{i}, 'corr');
            end
            test_array = CBIG_MMP_check_example_results(fullfile(outdir, 'results'), 1);
            
            % replace unit test results
            if replace_unit_test
                fprintf('Replacing unit test results... \n')
                copyfile(fullfile(outdir, 'results','KRR_corr_results.mat'), ...
                    fullfile(example_dir, 'ref_output'))
                copyfile(fullfile(outdir, 'results','LRR_corr_results.mat'), ...
                    fullfile(example_dir, 'ref_output'))
                copyfile(fullfile(outdir, 'results','Elasticnet_corr_results.mat'), ...
                    fullfile(example_dir, 'ref_output'))
                test_array = CBIG_MMP_check_example_results(fullfile(outdir, 'results'), 1);
            end
            rmpath(genpath(example_dir))

            assert(sum(test_array) == 6, 'Unit test failed: Examples do not match')
        end
        
        function test_example_second_level(testCase)
            % check whether second level models match: stacking and
            % multiKRR. First-level results must already be generated.
            
             % set directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            project_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Ooi2022_MMP');
            example_dir = fullfile(project_dir, 'examples');
            outdir = fullfile(project_dir, 'unit_tests', 'output');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            % run the example
            command = [example_dir, ...
                '/CBIG_MMP_HCP_example_multifeature_regression_wrapper.sh ' ...
                outdir];
            [status, cmdout] = system(command);

            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N Ooi2022_Example | wc -l'";
            [status, cmdout] = system(check_command);
            while str2num(cmdout) >= 1
                fprintf("Waiting for job(s) to finish... \n")
                pause(30)
                [status, cmdout] = system(check_command);
            end
            
            % compare the results
            addpath(genpath(example_dir))
            CBIG_MMP_HCP_collate_result_example_wrapper( ...
                fullfile(outdir, 'output'), fullfile(outdir, 'results'), ...
                'combined_models', 'corr');
            test_array = CBIG_MMP_check_example_results(fullfile(outdir, 'results'), 2);
            
            % replace unit test results
            if replace_unit_test
                fprintf('Replacing unit test results... \n')
                copyfile(fullfile(outdir, 'results','combined_models_corr_results.mat'), ...
                    fullfile(example_dir, 'ref_output'))
                test_array = CBIG_MMP_check_example_results(fullfile(outdir, 'results'), 1);
            end
            rmpath(genpath(example_dir))

            assert(sum(test_array) == 2, 'Unit test failed: Examples do not match')
            % remove the output directory
            rmdir(outdir, 's')
        end
    end
end

