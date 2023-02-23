classdef CBIG_ArealMSHBM_unit_test < matlab.unittest.TestCase
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2022_ArealMSHBM', 'unit_tests');
            out_dir = fullfile(cur_dir, 'output');
            mkdir(out_dir)
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2022_ArealMSHBM', 'examples'));
            
            CBIG_ArealMSHBM_example_wrapper(out_dir);
            
            % check example results
            CBIG_ArealMSHBM_check_example_results(out_dir, 1);
            
            % replace reference output if flag is 1
            ref_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation',...
                'Kong2022_ArealMSHBM', 'examples', 'ref_results');
            par = fullfile(out_dir, 'generate_individual_parcellations', 'ind_parcellation_gMSHBM', 'test_set',...
            '4_sess', 'beta5', 'Ind_parcellation_MSHBM_sub2_w50_MRF50_beta5.mat');
            par_ref = fullfile(ref_dir, 'generate_individual_parcellations', 'ind_parcellation_gMSHBM', 'test_set',...
            '4_sess', 'beta5', 'Ind_parcellation_MSHBM_sub2_w50_MRF50_beta5.mat');
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_MSHBM_unit_test...');
                copyfile(par, par_ref);
            end

            % remove the output directory
            rmdir(out_dir, 's')

            % remove path
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Kong2022_ArealMSHBM', 'examples'));
        end
    end
end