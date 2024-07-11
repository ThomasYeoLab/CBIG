classdef CBIG_MMM_unit_test < matlab.unittest.TestCase
% Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_KRR(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Chen2024_MMM');
            cur_dir = fullfile(base_dir, 'unit_tests');
            code_dir = fullfile(base_dir, 'KRR_CLASSICAL');

            input_dir = fullfile(cur_dir, 'data');
            out_dir = fullfile(cur_dir, 'output');
            output_classical_dir = fullfile(out_dir, 'output_KRR_classical_exp_test');
            subj_list = fullfile(input_dir, 'exp_test', 'exp_test_subj_list.txt');
            phe_csv = fullfile(input_dir, 'exp_test', 'exp_test_phe_tab.csv');
            FC_file = fullfile(input_dir, 'exp_test', 'exp_test_fc.mat');
            phe_list = fullfile(input_dir, 'exp_test', 'exp_test_phe_list.txt');
            temp = textscan(fopen(phe_list), '%s');
            phes = temp{1};
            ks = {'10', '20', '50', '100', '200'};
            rngs_num = 3;

            addpath(code_dir);
            for i = 1:length(phes)
                for j = 1:length(ks)
                    CBIG_MMM_KRR_classical(CBIG_CODE_DIR, subj_list, phe_csv, FC_file,...
                        output_classical_dir, rngs_num, phes{i}, ks{j}, false, 'exp_test');
                end
            end
            CBIG_MMM_KRR_classical_summary(CBIG_CODE_DIR, output_classical_dir, phe_list, rngs_num, 'exp_test');
            rmpath(code_dir);

            % get reference output folder
            ref_output_dir = fullfile(cur_dir, 'ref_results');
            addpath(genpath(fullfile(base_dir, 'examples')));
            % compare the results or replace the reference output
            load(fullfile(CBIG_CODE_DIR, 'unit_tests','replace_unittest_flag'));
            if replace_unittest_flag
                % replace regression results
                copyfile(fullfile(output_classical_dir, 'final_result', 'krr_classical_res.mat'),...
                    fullfile(ref_output_dir, 'output_KRR_classical_exp', 'final_result',...
                    'krr_classical_res.mat'));
            else
                CBIG_MMM_check_example_results(ref_output_dir, output_classical_dir);
            end
            rmpath(genpath(fullfile(base_dir, 'examples')));

            % remove the output directory and path
            rmdir(out_dir, 's')
        end
    end
end