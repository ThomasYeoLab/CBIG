classdef CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_ut_test < matlab.unittest.TestCase
    % Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    methods (Test)
        
        function Test(testCase)
            
            % main estimation function
            CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_ut();
            % load output
            cd 'unit_tests/step2_simulation/';
            load([pwd '/save/ref_output.mat']);
            load([pwd '/save/sim_output.mat']);
            
            % compare the current output with expected output
            assert(size(CC_save,1) == 1,'no. of rows must be 1')
            assert(size(CC_save,2) == 3,'no. of rows must be 3')
            assert(all(all(abs(CC_save-CC_check) < 10e-6)),'correlation result is not correct')
            
        end
        
    end    
    
end