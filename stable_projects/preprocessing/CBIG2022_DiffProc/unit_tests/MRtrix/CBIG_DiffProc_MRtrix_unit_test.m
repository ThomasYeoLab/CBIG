classdef CBIG_DiffProc_MRtrix_unit_test < matlab.unittest.TestCase
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    methods(Test)
        
        function test_iFOD(TestCase)
            % test MRtrix for probabilistic tractography
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'examples'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
                'CBIG2022_DiffProc','MRtrix_example_data');
            base_output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'MRtrix', 'output');
            output_dir = fullfile(base_output_dir, 'test_iFOD');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's');
            end
            mkdir(output_dir);
                       
            command = [fullfile(DIFF_CODE_DIR, 'MRtrix', 'CBIG_DiffProc_batch_tractography.sh'), ' -s ', ...
                fullfile(ref_dir, 'input', 'MRtrix_subs.txt'), ' -i ', ...
                fullfile(ref_dir, 'input'), ' -o ', fullfile(output_dir), ...
                ' -m ', fullfile(ref_dir, 'input', 'b0_mask'), ' -p CBIG_py3 -t 2M'];
            [status, cmdout] = system(command);
            
            % compare results
            subj_output_dir = fullfile(output_dir, 'iFOD2', 'output', ...
                'sub-NDARINV0A4P0LWM_ses-baselineYear1Arm1_run-01_dwi', 'connectomes');
            test_array = CBIG_DiffProc_MRtrix_check_example_results(subj_output_dir, 'test_iFOD2');
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')
                copyfile(fullfile(subj_output_dir), fullfile(ref_dir, 'ref_output', 'test_iFOD2'));
            end
            assert(sum(test_array) == 4, 'Unit test failed: Examples do not match')
            % remove output directory if unit test passes
            rmdir(output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'examples'));   
        end
        
        function test_FACT(TestCase)
            % test MRtrix for deterministic tractography
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'examples'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
                'CBIG2022_DiffProc','MRtrix_example_data');
            base_output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'MRtrix', 'output');
            output_dir = fullfile(base_output_dir, 'test_FACT');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's');
            end
            mkdir(output_dir);
            
            command = [fullfile(DIFF_CODE_DIR, 'MRtrix', 'CBIG_DiffProc_batch_tractography.sh'), ' -s ', ...
                fullfile(ref_dir, 'input', 'MRtrix_subs.txt'), ' -i ', ...
                fullfile(ref_dir, 'input'), ' -o ', fullfile(output_dir), ...
                ' -m ', fullfile(ref_dir, 'input', 'b0_mask'), ' -a FACT -p CBIG_py3 -t 2M'];
            [status, cmdout] = system(command);
            
            % compare results
            subj_output_dir = fullfile(output_dir, 'FACT', 'output', ...
                'sub-NDARINV0A4P0LWM_ses-baselineYear1Arm1_run-01_dwi', 'connectomes');
            test_array = CBIG_DiffProc_MRtrix_check_example_results(subj_output_dir, 'test_FACT');
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')
                copyfile(fullfile(subj_output_dir), fullfile(ref_dir, 'ref_output', 'test_FACT'));
            end
            assert(sum(test_array) == 4, 'Unit test failed: Examples do not match')
            % remove output directory if unit test passes
            rmdir(output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'examples'));  
            
            % remove all directories if final unit test passes
            output_master_dir = fullfile(base_output_dir);
            if numel(dir(output_master_dir)) == 2
                rmdir(output_master_dir, 's');
            end
        end
        
    end
    
end
