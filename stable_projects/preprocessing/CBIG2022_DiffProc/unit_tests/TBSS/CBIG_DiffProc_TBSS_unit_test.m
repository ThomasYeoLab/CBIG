classdef CBIG_DiffProc_TBSS_unit_test < matlab.unittest.TestCase
    % Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  
    methods(Test)
        
        function testFSLbasic(TestCase)
            % test the base FSL code 
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'examples'));
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
                'CBIG2022_DiffProc','TBSS');
            base_output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'TBSS', 'output');
            output_dir = fullfile(base_output_dir, 'fsl_non_parallel');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            mkdir(output_dir)
            
            % generate TBSS
            copyfile(fullfile(ref_dir,'input'), output_dir);
            
            cd(output_dir)
            command = 'tbss_1_preproc *.nii.gz';
            [status, cmdout] = system(command);
            command = 'tbss_2_reg -T';
            [status, cmdout] = system(command);
            command = 'tbss_3_postreg -S';
            [status, cmdout] = system(command);
            command = 'tbss_4_prestats 0.2';
            [status, cmdout] = system(command);
            command = 'tbss_non_FA MD';
            [status, cmdout] = system(command);
            cd(fullfile(DIFF_CODE_DIR, 'unit_tests'));
            
            % compare results
            test.FA_dir = fullfile(output_dir, 'stats', 'all_FA_skeletonised.nii.gz');
            test.MD_dir = fullfile(output_dir, 'stats', 'all_MD_skeletonised.nii.gz');
            test_array = CBIG_DiffProc_TBSS_check_example_results(test);
            
            % replace unit test if flag is 1
            if replace_unit_test     
                fprintf('Replacing reference results \n')
                copyfile(test.FA_dir, fullfile(ref_dir, 'ref_output', 'all_FA_skeletonised.nii.gz'));
                copyfile(test.MD_dir, fullfile(ref_dir, 'ref_output', 'all_MD_skeletonised.nii.gz'));
            end
            
            assert(sum(test_array) == 8, 'Unit test failed: Examples do not match')
            rmdir(output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'examples'));           
        end
        
        function testFSLparallel(TestCase)
            % test the parallelized version of FSL code 
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            CBIG_TEST_DIR = getenv('CBIG_TESTDATA_DIR');
            DIFF_CODE_DIR = fullfile(CBIG_CODE_DIR, 'stable_projects', 'preprocessing', ...
                'CBIG2022_DiffProc');
            addpath(fullfile(DIFF_CODE_DIR,'examples'));
            script_dir = fullfile(DIFF_CODE_DIR, 'TBSS');
            ref_dir = fullfile(CBIG_TEST_DIR,'stable_projects','preprocessing', ...
                'CBIG2022_DiffProc','TBSS');
            base_output_dir = fullfile(DIFF_CODE_DIR, 'unit_tests', 'TBSS', 'output');
            output_dir = fullfile(base_output_dir, 'fsl_parallel');
            replace_unit_test = load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', 'replace_unittest_flag'));
            if exist(output_dir, 'dir')
                rmdir(output_dir, 's')
            end
            mkdir(output_dir)
            
            % generate TBSS
            copyfile(fullfile(ref_dir,'input'), output_dir);
            % run wrapper           
            command = [fullfile(script_dir,'CBIG_DiffProc_TBSS_wrapper.sh '), ' ', output_dir];
            [status, cmdout] = system(command);
            
            % compare results
            test.FA_dir = fullfile(output_dir, 'stats', 'all_FA_skeletonised.nii.gz');
            test.MD_dir = fullfile(output_dir, 'stats', 'all_MD_skeletonised.nii.gz');
            test_array = CBIG_DiffProc_TBSS_check_example_results(test);
            
            % No replacement for unit test coded since reference is same
            % as basic unit test
            assert(sum(test_array) == 8, 'Unit test failed: Examples do not match')
            
            % remove all directories if final unit test passes
            rmdir(base_output_dir, 's');
            rmpath(fullfile(DIFF_CODE_DIR,'examples'));
        end
        
    end
    
end
