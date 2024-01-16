classdef CBIG_DiffProc_AMICO_unit_test < matlab.unittest.TestCase
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    methods(Test)
        
        function test_basic(TestCase)
            % test NODDI fitting using AMICO package
            
            % prepare directories
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'examples'));
            input_dir = fullfile(CBIG_TEST_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc', 'AMICO_example_data');
            ref_dir = fullfile(DIFF_CODE_DIR,'examples','AMICO_example_data');
            output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'AMICO', 'output');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            mkdir(output_dir)
            
            % generate AMICO
            command = [fullfile(DIFF_CODE_DIR, 'AMICO', 'CBIG_DiffProc_runAMICO.sh'), ' -s ', ...
                fullfile(input_dir, 'AMICO_subs.txt'), ...
                ' -d ', fullfile(input_dir, 'diffusion_data'), ' -o ', output_dir, ...
                ' -m ', fullfile(input_dir, 'b0_mask'), ' -p CBIG_py3'];
            [status, cmdout] = system(command);

            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N AMICO | wc -l'";
            [status, cmdout] = system(check_command);
            while str2num(cmdout) >= 1
                fprintf("Waiting for job to finish... \n")
                pause(30)
                [status, cmdout] = system(check_command);
            end
            
            % collate output paths
            subj_output_dir = fullfile(output_dir, 'output', 'AMICO', ...
                'subject_01', 'AMICO', 'NODDI');
            test.OD_dir = fullfile(subj_output_dir, 'FIT_OD.nii.gz');
            test.ICVF_dir = fullfile(subj_output_dir, 'FIT_ICVF.nii.gz');
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')   
                copyfile(fullfile(subj_output_dir, 'FIT_OD.nii.gz'), ...
                    fullfile(ref_dir, 'ref_output', 'output', 'AMICO', ...
                    'subject_01', 'AMICO', 'NODDI', 'FIT_OD.nii.gz'));
                copyfile(fullfile(subj_output_dir, 'FIT_ICVF.nii.gz'), ...
                    fullfile(ref_dir, 'ref_output', 'output', 'AMICO', ...
                    'subject_01', 'AMICO', 'NODDI', 'FIT_ICVF.nii.gz'));
            end
            % check if results are consistent
            test_array = CBIG_DiffProc_AMICO_check_example_results(test);
            assert(sum(test_array) == 4, 'Unit test failed: Examples do not match')
            % remove output directory if unit test passes
            rmdir(output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'examples'));   
        end
        
    end
    
end
