classdef CBIG_KRDNN_unit_test < matlab.unittest.TestCase
% Written by He Tong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_KR_TVT(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'unit_tests');
            % load replace_unittest_flag
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));

            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
            
            CBIG_KRDNN_KR_TVT_example();

            % get output folder
            out_dir = fullfile(cur_dir, 'output', 'output_kr_tvt_example');
            
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_KRDNN, test_KR_TVT...');
                ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                    'predict_phenotypes', 'He2019_KRDNN', 'examples', ...
                    'results', 'KR_TVT');
                ref_output = load(fullfile(ref_dir, 'final_result.mat'));
                test_output = load(fullfile(out_dir, 'final_result.mat'));
                % check coorelation
                abscorr = abs(ref_output.optimal_corr - test_output.optimal_corr);
                disp(['Total correlation difference (KRTVT) :' num2str(sum(sum(abscorr)))]);
                % check MAE of age prediction
                absmae = abs(ref_output.age_mae - test_output.age_mae);
                disp(['Total Age MAE differnce (KRTVT) :' num2str(sum(sum(absmae)))]);
                % check accuracy of sex prediction
                absacc = abs(ref_output.sex_accuracy - test_output.sex_accuracy);
                disp(['Total Sex prediction accuracy differnce (KRTVT) :' num2str(sum(sum(absacc)))]);
                ref_output = test_output;
                save(fullfile(ref_dir, 'final_result.mat'), '-struct', 'ref_output');
                
            else
                CBIG_KRDNN_check_example_results(out_dir, false);
            end
            % remove the output directory and path
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
        end

        function test_KR_CV(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'unit_tests');
            % load replace_unittest_flag
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
            
            CBIG_KRDNN_KR_CV_example();

            % get output folder
            out_dir = fullfile(cur_dir, 'output', 'output_kr_cv_example');
            
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_KRDNN, test_KR_TVT...');
                ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                    'predict_phenotypes', 'He2019_KRDNN', 'examples', ...
                    'results', 'KR_CV_matlab_r2018b');
                ref_output = load(fullfile(ref_dir, 'final_result.mat'));
                test_output = load(fullfile(out_dir, 'final_result.mat'));
                abserror = abs(ref_output.optimal_acc - test_output.optimal_acc);
                disp(['Total acc difference (KRCV) :' num2str(sum(sum(abserror)))]);
                ref_output = test_output;
                save(fullfile(ref_dir, 'final_result.mat'),  '-struct', 'ref_output');
                
            else
                CBIG_KRDNN_check_example_results(out_dir);
            end
            
            % remove the output directory and path
            rmdir(out_dir, 's')
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2019_KRDNN', 'examples'));
            rmdir(fullfile(cur_dir, 'output'))
        end
    end
end