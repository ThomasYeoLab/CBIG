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
            
            % compare the results or replace the reference output
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests','replace_unittest_flag'));
            if replace_unittest_flag
                % replace regression results
                copyfile(fullfile(outdir,'LRR','final_result_all_score.mat'), ...
                    fullfile(example_dir,'ref_output','LRR','final_result_all_score.mat'));
                copyfile(fullfile(outdir,'singleKRR','final_result_all_score.mat'), ...
                    fullfile(example_dir,'ref_output','singleKRR','final_result_all_score.mat'));
                copyfile(fullfile(outdir,'multiKRR','final_result_all_score.mat'), ...
                    fullfile(example_dir,'ref_output','multiKRR','final_result_all_score.mat'));
                % replace PFM results
                copyfile(fullfile(outdir,'PFM','LRR','PFM_score1_all_folds.mat'), ...
                    fullfile(example_dir,'ref_output','PFM','LRR','PFM_score1_all_folds.mat'));
                copyfile(fullfile(outdir,'PFM','LRR','PFM_score2_all_folds.mat'), ...
                    fullfile(example_dir,'ref_output','PFM','LRR','PFM_score2_all_folds.mat'));

                copyfile(fullfile(outdir,'PFM','singleKRR','PFM_score1_all_folds.mat'), ...
                    fullfile(example_dir,'ref_output','PFM','singleKRR','PFM_score1_all_folds.mat'));
                copyfile(fullfile(outdir,'PFM','singleKRR','PFM_score2_all_folds.mat'), ...
                    fullfile(example_dir,'ref_output','PFM','singleKRR','PFM_score2_all_folds.mat'));

                copyfile(fullfile(outdir,'PFM','multiKRR','PFM_score1_all_folds.mat'), ...
                    fullfile(example_dir,'ref_output','PFM','multiKRR','PFM_score1_all_folds.mat'));
                copyfile(fullfile(outdir,'PFM','multiKRR','PFM_score2_all_folds.mat'), ...
                    fullfile(example_dir,'ref_output','PFM','multiKRR','PFM_score2_all_folds.mat'));
            else
                CBIG_TRBPC_check_example_results(outdir)
            end

            % remove the output directory
            rmdir(outdir, 's')
            rmpath(example_dir)
        end
    end
end