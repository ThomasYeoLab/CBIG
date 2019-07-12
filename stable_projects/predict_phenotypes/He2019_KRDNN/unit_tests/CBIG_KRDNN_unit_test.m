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
            
            % check the results
            ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples', ...
                'results', 'KR_TVT');

            % get output folder
            out_dir = fullfile(cur_dir, 'output_kr_tvt_example');
            
            % % compare prediction
            ref_output = load(fullfile(ref_dir, 'final_result.mat'));
            test_output = load(fullfile(out_dir, 'final_result.mat'));

            % check correlation
            diff_corr = sum(abs(ref_output.optimal_corr - test_output.optimal_corr));
            assert(diff_corr < 1e-6, sprintf('optimal_corr difference: %f', diff_corr));

            % check MAE of age prediction
            diff_mae = abs(ref_output.age_mae - test_output.age_mae);
            assert(diff_mae < 1e-6, sprintf('age_mae difference: %f', diff_mae));

            % check accuracy of sex prediction
            diff_acc = abs(ref_output.sex_accuracy - test_output.sex_accuracy);
            assert(diff_acc < 1e-6, sprintf('sex_accuracy difference: %f', diff_acc));
            
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
            
            % check the results
            ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples', ...
                'results', 'KR_CV_matlab_r2018b');

            % get output folder
            out_dir = fullfile(cur_dir, 'output_kr_cv_example');
            
            % % compare prediction
            ref_output = load(fullfile(ref_dir, 'final_result.mat'));
            test_output = load(fullfile(out_dir, 'final_result.mat'));

            % check accuracy
            diff_acc = sum(sum(abs(ref_output.optimal_acc - test_output.optimal_acc)));
            assert(diff_acc < 1e-6, sprintf('optimal_acc difference: %f', diff_acc));
            
            % remove the output directory and path
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
        end
    end
end