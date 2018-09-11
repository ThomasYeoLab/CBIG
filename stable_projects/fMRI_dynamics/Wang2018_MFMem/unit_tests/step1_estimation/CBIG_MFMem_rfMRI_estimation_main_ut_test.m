classdef CBIG_MFMem_rfMRI_estimation_main_ut_test < matlab.unittest.TestCase
    % Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
    methods (Test)
        
        function Test(testCase)
            
            % main estimation function
            CBIG_MFMem_rfMRI_estimation_main_ut(30,3);
            % load output
            cd 'unit_tests/step1_estimation/';
            load([pwd '/save/corr_check.mat']);
            load([pwd '/save/corr_saved.mat']);
            
            % compare the current output with expected output
            assert(size(rrr_z,1) == 1,'no. of rows must be 1')
            assert(size(rrr_z,2) == 3,'no. of rows must be 3')
            assert(all(all(abs(rrr_z-rrr_z_check) < 10e-6)),'correlation result is not correct')
            
        end
        
    end    
    
end