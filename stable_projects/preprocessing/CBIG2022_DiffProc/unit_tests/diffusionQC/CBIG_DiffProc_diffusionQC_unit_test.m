classdef CBIG_DiffProc_diffusionQC_unit_test < matlab.unittest.TestCase
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    methods(Test)
        
        function test_diffusionQC_basic(TestCase)
            % basic test case
            % prepare directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); 
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc', 'diffusionQC');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'diffusionQC', 'output', 'basic');
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            
            % run QC pipeline
            command = strcat(fullfile(DIFF_CODE_DIR, 'CBIG_DiffProc_diffusionQC.sh'), ...
                " -s sub-1 -d ", fullfile(ref_dir, 'input', 'sub-1', 'ses-baselineYear1Arm1', 'dwi'), ...
                " -i sub-1_ses-baselineYear1Arm1_run-01_dwi.nii", ...
                " -a sub-1_ses-baselineYear1Arm1_run-01_dwi.bval -e sub-1_ses-baselineYear1Arm1_run-01_dwi.bvec ", ...
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
                  
            % read QC file and compare to reference
            subj_qc = fullfile(output_dir, 'QC_output', 'sub-1_QC.txt');
            fid = fopen(subj_qc,'r'); % variable names text file
            tline = fgetl(fid);
            i = 1;
            while ischar(tline)
                tline = fgetl(fid);
                try
                    test_qc_all{:,i} = split(tline, ": ");
                    test.qc_vals(:,i) = str2num(test_qc_all{:,i}{2});
                end
                i = i + 1;
            end
            fclose(fid);
            
            test.FA_dir = fullfile(output_dir, 'fdt', 'dti_FA.nii.gz');
            test.SSE_dir = fullfile(output_dir, 'fdt', 'dti_sse.nii.gz')
            test_array = CBIG_DiffProc_diffusionQC_check_unit_test_results(test, "basic");
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')
                qc_vals = test.qc_vals;
                save(fullfile(ref_dir, 'ref_output', 'diffusionQC_sub-1_basic.mat'),'qc_vals')
                copyfile(test.FA_dir, fullfile(ref_dir, 'ref_output', 'dti_FA.nii.gz'));
                copyfile(test.SSE_dir, fullfile(ref_dir, 'ref_output', 'dti_sse.nii.gz'));
            end
            assert(sum(test_array) == 5, 'Unit test failed: Examples do not match')
            rmdir(output_dir, 's');
            
            rmpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
        end
    
        function test_diffusionQC_rounding(TestCase)
            % test case where b values need to be rounded off to nearest
            % 1000
            % prepare directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); 
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc', 'diffusionQC');
            output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'diffusionQC', 'output', 'rounding');
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            % run QC pipeline
            command = strcat(fullfile(DIFF_CODE_DIR, 'CBIG_DiffProc_diffusionQC.sh'), " -s sub-1 -d ", ...
                fullfile(ref_dir, 'input', 'sub-1', 'ses-baselineYear1Arm1', 'dwi'), ...
                " -i sub-1_ses-baselineYear1Arm1_run-01_dwi.nii", ...
                " -a sub-1_ses-baselineYear1Arm1_run-01_dwi_mod.bval", ...
                " -e sub-1_ses-baselineYear1Arm1_run-01_dwi.bvec ", ...
                " -o ", output_dir, " -r");
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
                  
            % read QC file and compare to reference
            subj_qc = fullfile(output_dir, 'QC_output', 'sub-1_QC.txt');
            fid = fopen(subj_qc,'r'); % variable names text file
            tline = fgetl(fid);
            i = 1;
            while ischar(tline)
                tline = fgetl(fid);
                try
                    test_qc_all{:,i} = split(tline, ": ");
                    test.qc_vals(:,i) = str2num(test_qc_all{:,i}{2});
                end
                i = i + 1;
            end
            fclose(fid);
            
            
            test_array = CBIG_DiffProc_diffusionQC_check_unit_test_results(test, "rounding");
            % replace unit test if flag is 1
            if replace_unit_test
                fprintf('Replacing reference results \n')
                qc_vals = test.qc_vals;
                save(fullfile(ref_dir, 'ref_output', 'diffusionQC_sub-1_round.mat'),'qc_vals')
            end
            assert(sum(test_array) == 1, 'Unit test failed: Examples do not match')
            rmdir(output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
        end
    
        function test_diffusionQC_single_shell(TestCase)
            % test case where only b=1000 shell is used
            % prepare directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); 
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc', 'diffusionQC');
            output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'diffusionQC', 'output', 'single_shell');
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            
            % run QC pipeline
            command = strcat(fullfile(DIFF_CODE_DIR, 'CBIG_DiffProc_diffusionQC.sh'), ...
                " -s sub-1 -d ", fullfile(ref_dir, 'input', 'sub-1', 'ses-baselineYear1Arm1', 'dwi'), ...
                " -i sub-1_ses-baselineYear1Arm1_run-01_dwi.nii", ...
                " -a sub-1_ses-baselineYear1Arm1_run-01_dwi.bval -e sub-1_ses-baselineYear1Arm1_run-01_dwi.bvec ", ...
                " -o ", output_dir, " -u 1000");
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); 
            submit_command = fullfile(CBIG_CODE_DIR,'setup', strcat("CBIG_pbsubmit -cmd '", ...
                command, "' -walltime '00:30:00' -name 'DiffQC_unittest'"));
            [status, cmdout] = system(submit_command);
            
            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N DiffQC_unittest | wc -l'";
            [status, cmdout] = system(check_command);
            while str2num(cmdout) >= 1
                fprintf("Waiting for job to finish... \n")
                pause(30)
                [status, cmdout] = system(check_command);
            end
                  
            % read QC file and compare to reference
            subj_qc = fullfile(output_dir, 'QC_output', 'sub-1_QC.txt');
            fid = fopen(subj_qc,'r'); % variable names text file
            tline = fgetl(fid);
            i = 1;
            while ischar(tline)
                tline = fgetl(fid);
                try
                    test_qc_all{:,i} = split(tline, ": ");
                    test.qc_vals(:,i) = str2num(test_qc_all{:,i}{2});
                end
                i = i + 1;
            end
            fclose(fid);
            
            test.bval_dir = fullfile(output_dir, 'fdt', 'sub-1_single_shell_1000.bval');
            test.bvec_dir = fullfile(output_dir, 'fdt', 'sub-1_single_shell_1000.bvec');
            test.single_shell_dir = fullfile(output_dir, 'fdt', 'sub-1_single_shell_1000.nii.gz');
            test.FA_dir = fullfile(output_dir, 'fdt', 'dti_FA.nii.gz');
            test.SSE_dir = fullfile(output_dir, 'fdt', 'dti_sse.nii.gz')
            test_array = CBIG_DiffProc_diffusionQC_check_unit_test_results(test, "single_shell");
            
            % replace unit test if flag is 1
            if replace_unit_test    
                fprintf('Replacing reference results \n')
                qc_vals = test.qc_vals;
                save(fullfile(ref_dir, 'ref_output', 'diffusionQC_sub-1_single_shell.mat'),'qc_vals')
                copyfile(test.FA_dir, fullfile(ref_dir, 'ref_output', 'dti_FA_ss.nii.gz'));
                copyfile(test.SSE_dir, fullfile(ref_dir, 'ref_output', 'dti_sse_ss.nii.gz'));
                copyfile(test.bval_dir, fullfile(ref_dir, 'ref_output', 'sub-1_single_shell_1000.bval'));
                copyfile(test.bvec_dir, fullfile(ref_dir, 'ref_output', 'sub-1_single_shell_1000.bvec'));
                copyfile(test.single_shell_dir, ...
                    fullfile(ref_dir, 'ref_output', 'sub-1_single_shell_1000.nii.gz'));
            end
            assert(sum(test_array) == 9, 'Unit test failed: Examples do not match')
            rmdir(output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'unit_tests'));
             
            % remove all directories if final unit test passes
            output_master_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'diffusionQC', 'output');
            if numel(dir(output_master_dir)) == 2
                rmdir(output_master_dir, 's');
            end
        end
    end
end
