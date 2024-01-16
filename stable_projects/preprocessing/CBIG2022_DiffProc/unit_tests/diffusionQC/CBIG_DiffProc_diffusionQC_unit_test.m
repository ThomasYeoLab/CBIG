classdef CBIG_DiffProc_diffusionQC_unit_test < matlab.unittest.TestCase
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    methods(Test)
        
        function test_diffusionQC_basic(TestCase)
            % basic test case: round b values to nearest 100 and generate
            % diffusion tensor results
            
            % prepare directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); 
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'examples'));
            input_dir = fullfile(CBIG_TEST_DIR,'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc', 'diffusionQC_example_data');
            ref_dir = fullfile(DIFF_CODE_DIR,'examples','diffusionQC_example_data');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'diffusionQC', 'output');
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            
            % run QC pipeline
            command = strcat(fullfile(DIFF_CODE_DIR, 'CBIG_DiffProc_diffusionQC.sh'), ...
                " -s sub-A00059845 -d ", fullfile(input_dir, 'input', 'sub-A00059845', 'ses-DS2', 'dwi'), ...
                " -i sub-A00059845_ses-DS2_dwi.nii.gz -a sub-A00059845_ses-DS2_dwi.bval", ...
                " -e sub-A00059845_ses-DS2_dwi.bvec --round_bvals_100 ", ...
                " -o ", output_dir);
            submit_command = strcat(fullfile(CBIG_CODE_DIR, 'setup', 'CBIG_pbsubmit'), " -cmd '", ...
                command, "' -walltime '0:30:00' -name 'DiffQC_unittest' -mem '8G'");
            [status, cmdout] = system(submit_command);

            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N DiffQC_unittest | wc -l'";
            [status, cmdout] = system(check_command);
            while str2num(cmdout) >= 1
                fprintf("Waiting for job to finish... \n")
                pause(30)
                [status, cmdout] = system(check_command);
            end
            
            % collate output paths
            test.qc_txt = fullfile(output_dir, 'QC_output', 'sub-A00059845_QC.txt');
            test.FA_dir = fullfile(output_dir, 'fdt', 'dti_FA.nii.gz');
            test.SSE_dir = fullfile(output_dir, 'fdt', 'dti_sse.nii.gz')
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')
                copyfile(test.qc_txt, fullfile(ref_dir,'ref_output','QC_output','sub-A00059845_QC.txt'));
                copyfile(test.FA_dir, fullfile(ref_dir,'ref_output','fdt','dti_FA.nii.gz'));
                copyfile(test.SSE_dir, fullfile(ref_dir,'ref_output','fdt','dti_sse.nii.gz'));
            end
            % check if results are consistent
            test_array = CBIG_DiffProc_diffusionQC_check_example_results(test);
            assert(sum(test_array) == 5, 'Unit test failed: Examples do not match')
            rmdir(output_dir, 's');
            
            rmpath(fullfile(DIFF_CODE_DIR,'examples'));
        end
        
    end
end
