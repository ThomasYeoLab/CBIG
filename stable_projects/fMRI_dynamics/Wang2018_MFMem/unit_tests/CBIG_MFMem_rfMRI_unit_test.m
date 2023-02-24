classdef CBIG_MFMem_rfMRI_unit_test < matlab.unittest.TestCase
    % In this script, we perform the step_1_estimation unit test of Wang2018_MFMem project
    
    % The content of this script is bacically the same as this file's: CBIG_MFMem_rfMRI_estimation_main_ut_test
    % Here we change the test scrip's name and put it under our repo's unit_tests folder, 
    % so that it's easier for our wrapper script to call it
    
    % Written by Xiaolu Kong, Yang-Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        
        function Test_Estimation(testCase)            
            %% path setting
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'fMRI_dynamics', 'Wang2018_MFMem', 'unit_tests');
            ReferenceDir = fullfile(UnitTestDir, 'step1_estimation/save');
            OutputDir = fullfile(UnitTestDir, 'output');
            
            % create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call CBIG_MFMem_rfMRI_estimation_main_ut to generate results
            % add path
            addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'fMRI_dynamics', 'Wang2018_MFMem')));
            
            % main estimation function
            CBIG_MFMem_rfMRI_estimation_main_ut(30, 3, OutputDir);
            
            
            %% compare the current output with expected output
            load(fullfile(getenv('CBIG_CODE_DIR'), 'unit_tests', ...
                'replace_unittest_flag'))
            
            % load output
            load(fullfile(ReferenceDir, 'corr_check.mat')); % expected output
            load(fullfile(OutputDir, 'corr_saved.mat')); % new output
            
            if(replace_unittest_flag)
                disp('Replacing unit test reference results of MFMem, estimation case...')
                rrr_z_check = rrr_z;
                save(fullfile(ReferenceDir, 'corr_check.mat'),'rrr_z_check')
            else
                assert(size(rrr_z,1) == 1, 'no. of rows must be 1')
                assert(size(rrr_z,2) == 3, 'no. of rows must be 3')
                assert(all(all(abs(rrr_z-rrr_z_check) < 10e-6)), ...
                    'correlation result is not correct')
            end
            
            %% remove path and intermediate output data (IMPORTANT)
            rmpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
                'fMRI_dynamics', 'Wang2018_MFMem')));
            rmdir(OutputDir, 's');
            
        end
        
    end    
    
end