classdef CBIG_pFIC_unit_test < matlab.unittest.TestCase
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            load(fullfile(CBIG_CODE_DIR, 'unit_tests', 'replace_unittest_flag'));
            output_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'fMRI_dynamics', 'Zhang2024_pFIC', 'examples');
            output_dir = fullfile(output_dir, 'output');
            mkdir(output_dir)
            
            % run the example
            system('./CBIG_pFIC_unit_test_submission_matlab.sh')
            
            % check if job finished
            check_command = "ssh headnode 'qselect -u $(whoami) -N pFIC_unittest | wc -l'";
            [status, cmdout] = system(check_command);
            timer = 0;
            while str2num(cmdout) >= 1
                fprintf("Waiting for the job to finish... \n")
                timer = timer + 150;
                if timer > 3600
                    disp('Taking too long to finish... Consider running it manually.')
                    break
                end                    
                pause(150)
                [status, cmdout] = system(check_command);
            end
            
            % check example results
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'fMRI_dynamics', 'Zhang2024_pFIC', 'examples'));            
            
            reference_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'fMRI_dynamics', 'Zhang2024_pFIC', 'examples', 'reference_output');
            CBIG_pFIC_check_example_results(output_dir, reference_dir);
            % replace reference output if flag is 1
            if(replace_unittest_flag)
                disp('Replacing unit test reference results for CBIG_pFIC_unit_test...');
                copyfile(output_dir, reference_dir);
            end

            % remove the output directory
            rmdir(output_dir, 's')

            % remove path
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'fMRI_dynamics', 'Zhang2024_pFIC', 'examples'));
        end
    end
end
