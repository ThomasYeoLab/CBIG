classdef CBIG_MM_unit_test < matlab.unittest.TestCase
% Written by He Tong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_MM_KRR(testCase)
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'predict_phenotypes', 'He2022_MM');
            cur_dir = fullfile(base_dir, 'unit_tests');
            code_dir = fullfile(base_dir, 'KRR_CLASSICAL');

            input_dir = fullfile(cur_dir, 'data');
            out_dir = fullfile(cur_dir, 'output');
            output_classical_dir = fullfile(out_dir, 'output_KRR_classical_exp');
            subj_list = fullfile(input_dir, 'exp_test_subj_list.txt');
            phe_csv = fullfile(input_dir, 'exp_test_final.csv');
            FC_file = fullfile(input_dir, 'exp_test_pfc.mat');
            phe_list = fullfile(input_dir, 'exp_test_final_phe_list.txt');
            temp = textscan(fopen(phe_list), '%s');
            phes = temp{1};
            ks = {'10', '20', '50', '100', '200'};
            rngs_num = 2;
            
            addpath(code_dir);
            for i = 1:length(phes)
                for j = 1:length(ks)
                    CBIG_MM_KRR_classical(CBIG_CODE_DIR, subj_list, phe_csv, FC_file,...
                        output_classical_dir, rngs_num, phes{i}, ks{j}, false);
                end
            end
            CBIG_MM_KRR_classical_summary(CBIG_CODE_DIR, output_classical_dir, phe_list, rngs_num, 'exp');
            rmpath(code_dir);

            code_dir = fullfile(base_dir, 'KRR_MM');
            output_mm_dir = fullfile(out_dir, 'output_KRR_mm');
            phe_list = fullfile(input_dir, 'exp_train_final_phe_list.txt');
            phe_csv = fullfile(input_dir, 'exp_train_test_final.csv');
            subj_list = fullfile(input_dir, 'exp_train_test_subj_list.txt');
            fc_mat = fullfile(input_dir, 'exp_train_test_pfc.mat');
            subj_list_train = fullfile(input_dir, 'exp_train_subj_list.txt');
            subj_list_extra = fullfile(input_dir, 'exp_test_subj_list.txt');
            rngs_num = 1; % run the KRR MM base 1 times

            addpath(code_dir);
            for rng_num = 1:rngs_num
                CBIG_MM_KRR_MM_base(CBIG_CODE_DIR, subj_list, phe_csv, fc_mat, output_mm_dir,...
                    num2str(rng_num), phe_list, subj_list_train, subj_list_extra);
            end            
            mm_rng_nums = 2;
            prefix = 'exp';
            CBIG_MM_KRR_MM_summary(CBIG_CODE_DIR, base_dir, output_mm_dir, output_classical_dir,...
                input_dir, mm_rng_nums, prefix);
            rmpath(code_dir);

            % get reference output folder
            ref_output_dir = fullfile(cur_dir, 'results');
            addpath(genpath(fullfile(base_dir, 'examples')));
            % compare the results or replace the reference output
            load(fullfile(CBIG_CODE_DIR, 'unit_tests','replace_unittest_flag'));
            if replace_unittest_flag
                % replace regression results
                copyfile(fullfile(output_classical_dir, 'final_result', 'krr_classical_res_test.mat'),...
                    fullfile(ref_output_dir, 'output_KRR_classical_exp', 'final_result',...
                    'krr_classical_res_test.mat'));
                copyfile(fullfile(output_mm_dir, 'meta_result', 'meta_res_test.mat'), ...
                    fullfile(ref_output_dir, 'output_KRR_mm', 'meta_result', 'meta_res_test.mat'));
            else
                CBIG_MM_check_example_results(ref_output_dir, output_classical_dir, output_mm_dir);
            end
            rmpath(genpath(fullfile(base_dir, 'examples')));
            
            % remove the output directory and path
            rmdir(out_dir, 's')
        end
    end
end