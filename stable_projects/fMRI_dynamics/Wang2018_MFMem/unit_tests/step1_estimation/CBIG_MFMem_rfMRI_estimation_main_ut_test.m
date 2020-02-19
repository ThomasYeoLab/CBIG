% classdef CBIG_MFMem_rfMRI_estimation_main_ut_test < matlab.unittest.TestCase
% % Written by Xiaolu Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%     
% 
%     methods (Test)
%         
%         function Test_Estimation(testCase)            
%             %% path setting
%             UnitTestDir = [getenv('CBIG_CODE_DIR') '/unit_tests'];
%             FolderStructure = '/stable_projects/fMRI_dynamics/Wang2018_MFMem/';
%             ReferenceDir = [getenv('CBIG_CODE_DIR'), '/stable_projects/fMRI_dynamics/Wang2018_MFMem',...
%                 '/unit_tests/step1_estimation/save'];
%             OutputDir = fullfile(ReferenceDir, 'output'); % this output dir is case specific
%             
%             % create output dir (IMPORTANT)
%             if(exist(OutputDir, 'dir'))
%                 rmdir(OutputDir, 's')
%             end
%             mkdir(OutputDir);
%             
%             
%             %% call CBIG_MFMem_rfMRI_estimation_main_ut to generate results
%             % add path
%             addpath(genpath([getenv('CBIG_CODE_DIR'), '/stable_projects/fMRI_dynamics/Wang2018_MFMem']));
%             
%             % main estimation function
%             CBIG_MFMem_rfMRI_estimation_main_ut(30, 3, OutputDir);
%             
%             
%             %% compare the current output with expected output
%             % load output
%             load([ReferenceDir, '/corr_check.mat'],'rrr_z_check'); % expected output
%             load([OutputDir, '/corr_saved.mat'],'rrr_z'); % new output
%             
%             assert(size(rrr_z,1) == 1, 'no. of rows must be 1')
%             assert(size(rrr_z,2) == 3, 'no. of rows must be 3')
%             assert(all(all(abs(rrr_z-rrr_z_check) < 10e-6)), 'correlation result is not correct')
%             
%             
%             %% remove path and intermediate output data (IMPORTANT)
%             rmdir(OutputDir, 's');
%             rmpath(genpath([getenv('CBIG_CODE_DIR'), '/stable_projects/fMRI_dynamics/Wang2018_MFMem']));
%         end
%         
%     end    
%     
% end