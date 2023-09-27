classdef CBIG_DiffProc_AMICO_unit_test < matlab.unittest.TestCase
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    methods(Test)
        
        function test_singleslice(TestCase)
            % test NODDI fitting for a single slice
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
                'CBIG2022_DiffProc','AMICO');
            base_output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'AMICO', 'output');
            output_dir = fullfile(base_output_dir, 'test_singleslice');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            mkdir(output_dir)
                       
            command = [fullfile(DIFF_CODE_DIR, 'AMICO', 'CBIG_DiffProc_runAMICO.sh'), ' -s ', ...
                fullfile(ref_dir, 'input', 'AMICO_subs.txt'), ...
                ' -d ', fullfile(ref_dir, 'input'), ' -o ', fullfile(output_dir), ...
                ' -m ', fullfile(ref_dir, 'input', 'b0_mask'), ' -p CBIG_py3'];
            [status, cmdout] = system(command);

            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N AMICO | wc -l'";
            [status, cmdout] = system(check_command);
            while str2num(cmdout) >= 1
                fprintf("Waiting for job to finish... \n")
                pause(30)
                [status, cmdout] = system(check_command);
            end
            
            % compare results
            subj_output_dir = fullfile(output_dir, 'output', 'AMICO', ...
                'sub-1_ses-baselineYear1Arm1_run-01_dwi', 'AMICO', 'NODDI');
            test.OD_dir = fullfile(subj_output_dir, 'FIT_OD.nii.gz');
            test.ICVF_dir = fullfile(subj_output_dir, 'FIT_ICVF.nii.gz');
            % check example script prints if there are any differences
            test_array = CBIG_DiffProc_AMICO_check_unit_test_results(test);
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')              
                copyfile(fullfile(subj_output_dir, 'FIT_OD.nii.gz'), ...
                    fullfile(ref_dir, 'ref_output', 'FIT_OD.nii.gz'));
                copyfile(fullfile(subj_output_dir, 'FIT_ICVF.nii.gz'), ...
                    fullfile(ref_dir, 'ref_output', 'FIT_ICVF.nii.gz'));
            end
            assert(sum(test_array) == 4, 'Unit test failed: Examples do not match')
            % remove output directory if unit test passes
            rmdir(base_output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'unit_tests'));   
        end
        
    end
    
end
