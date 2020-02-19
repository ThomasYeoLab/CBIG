classdef CBIG_gwMRF_unit_test < matlab.unittest.TestCase
%
% Target project:
%                 Schaefer2018_LocalGlobal
%
% Case design:
%                 Schaefer2018_LocalGlobal already has unit test script, here 
%                 we call its unit test script and automatically
%                 make judgement on whether the unit test is passed or
%                 failed
%
%                 For now, the stable projects' unit tests in our repo
%                 require manual check of the output txt files or images,
%                 making it incovenient for wrapper function to call them
%                 and automatically draw conclusions. As a result, we write
%                 some simple matlab test functions for these stable
%                 projects' unit tests.
%
% Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_Case(testCase)
            %% path setting
            UnitTestDir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
                'brain_parcellation', 'Schaefer2018_LocalGlobal', 'unit_tests');
            OutputDir = fullfile(UnitTestDir, 'output'); 
            
            % create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);
            
            
            %% call Schaefer2018_LocalGlobal unit test script to generate results
            cmd = ['sh ' fullfile(UnitTestDir, 'scripts', 'CBIG_gwMRF_unit_test.sh'),...
                ' ', OutputDir];            
            system(cmd); % this will submit a job to HPC
            
            
            %% check output log file
            logfile_path = fullfile(OutputDir, 'logs', 'CBIG_gwMRF_unit_test.log');
            [~, error_messages] = system(['cat ', logfile_path, ' | grep FAILED']);
            [~, success_messages] = system(['cat ', logfile_path, ' | grep SUCCESS']);
            
            % extract detailed error meassages from log file
            assert(isempty(error_messages), sprintf(error_messages));
            assert(~isempty(success_messages), sprintf(success_messages));            
            
            %% remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
        end
        
        
    end
end