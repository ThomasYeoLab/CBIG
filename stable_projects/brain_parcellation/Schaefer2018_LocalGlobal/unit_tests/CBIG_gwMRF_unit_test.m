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
            
            %% prep work
            % create output dir (IMPORTANT)
            if(exist(OutputDir, 'dir'))
                rmdir(OutputDir, 's')
            end
            mkdir(OutputDir);

            %% generate example results
            % create example outdir, which is under the unit test parent folder
            ExampleDir = fullfile(OutputDir,'example_out');
            mkdir(ExampleDir);
            % generate example input file
            cmd = ['sh ${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'...
            '/examples/example_input/CBIG_gwMRF_create_example_input_fullpaths.sh  ' ExampleDir];
            system(cmd);
            
            %% run the rest of the unit test that checks intermediate results
            % generate unit test input full paths
            cmd = ['sh ' fullfile(UnitTestDir, 'scripts', 'CBIG_gwMRF_create_unit_tests_input_fullpaths.sh'),...
            ' ', OutputDir];            
            system(cmd);
            % call Schaefer2018_LocalGlobal unit test script to generate results (including example)
            % as well as check the result correctness. Its return value indicates whether unit test passed or failed.
            cmd = ['sh ' fullfile(UnitTestDir, 'scripts', 'CBIG_gwMRF_unit_test.sh'),...
                ' ', OutputDir];
            status = system(cmd); 
            assert(isequal(status,0), 'Unit test failed.');

            % check output for the example. This has to be done after we run CBIG_gwMRF_unit_test.sh,
            % since generating example results is done by this shell script. We choose to do this to
            % avoid switching between different environments for multiple times.
            addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','brain_parcellation',...
            'Schaefer2018_LocalGlobal','examples','scripts'));
            assert(CBIG_gwMRF_check_example_results(ExampleDir, 1),...
                'The example result for seed 1 does not match the reference result.');
            rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','brain_parcellation',...
            'Schaefer2018_LocalGlobal','examples','scripts'));
            
            %% remove intermediate output data (IMPORTANT)
            rmdir(OutputDir, 's');
            close(gcf);
        end  
    end
end