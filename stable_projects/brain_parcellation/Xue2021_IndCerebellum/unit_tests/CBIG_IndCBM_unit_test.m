classdef CBIG_IndCBM_unit_test < matlab.unittest.TestCase
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Xue2021_IndCerebellum', 'unit_tests');
            output_dir = fullfile(cur_dir, 'output');
            
            % run the example
            example_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'brain_parcellation', 'Xue2021_IndCerebellum', 'examples');
            addpath(example_dir);
            
            CBIG_IndCBM_example_wrapper(output_dir);
            
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_IndCBM_unit_test...');
                CBIG_IndCBM_check_example_results(output_dir);
                out_file1 = fullfile(output_dir, 'parcellation', 'sub1', 'IndCBM_parcellation_top100.dlabel.nii');
                out_file2 = fullfile(output_dir, 'parcellation', 'sub2', 'IndCBM_parcellation_top100.dlabel.nii');
                ref_file1 = fullfile(example_dir, 'ref_output', 'sub1', 'IndCBM_parcellation_top100.dlabel.nii');
                ref_file2 = fullfile(example_dir, 'ref_output', 'sub2', 'IndCBM_parcellation_top100.dlabel.nii');
                copyfile(out_file1, ref_file1);
                copyfile(out_file2, ref_file2);
            else
                % check example results
                assert(CBIG_IndCBM_check_example_results(output_dir), sprintf('Result check failed.'));
            end
                
            % remove the output directory
            rmdir(output_dir, 's')

            % remove path
            rmpath(example_dir);
        end
    end
end
