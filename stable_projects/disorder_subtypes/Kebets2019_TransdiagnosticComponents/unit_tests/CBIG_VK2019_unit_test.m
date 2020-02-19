classdef CBIG_VK2019_unit_test < matlab.unittest.TestCase
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            
            % create output folder
            current_dir = fileparts(mfilename('fullpath'));
            out_dir = fullfile(current_dir,'output');
            mkdir(out_dir)
            
            % run the example
            pos_v = strfind(current_dir,filesep);
            root_dir = fullfile(current_dir(1:pos_v(length(pos_v)) - 1));
            
            addpath(fullfile(root_dir,'examples'));
            
            CBIG_VK2019_example_wrapper(out_dir);
            
            % check the results
            ref_dir = fullfile(root_dir,'examples','correct_output');
            
            % compare RSFC & behavioral scores & loadings
            ref_Params = load(fullfile(ref_dir,'PLSresults_example.mat'));          
            ref_Lx = ref_Params.Lx;
            ref_Ly = ref_Params.Ly;
            ref_RSFC_loadings = ref_Params.LC_RSFC_loadings;
            ref_behav_loadings = ref_Params.LC_behav_loadings;
            
            test_Params = load(fullfile(out_dir,'PLSresults_example.mat'));
            test_Lx = test_Params.Lx;
            test_Ly = test_Params.Ly;
            test_RSFC_loadings = test_Params.LC_RSFC_loadings;
            test_behav_loadings = test_Params.LC_behav_loadings;
            
            diff_RSFC_scores = max(max(ref_Lx - test_Lx));
            diff_behav_scores = max(max(ref_Ly - test_Ly));
            diff_RSFC_loadings = max(max(ref_RSFC_loadings - test_RSFC_loadings));
            diff_behav_loadings = max(max(ref_behav_loadings - test_behav_loadings));
                        
            assert(diff_RSFC_scores < 1e-6, sprintf('maximum difference in RSFC scores: %f',diff_RSFC_scores));
            assert(diff_behav_scores < 1e-6, sprintf('maximum difference in behav scores: %f',diff_behav_scores));
            assert(diff_RSFC_loadings < 1e-6, sprintf('maximum difference in RSFC loadings: %f',diff_RSFC_loadings));
            assert(diff_behav_loadings < 1e-6, sprintf('maximum difference in behav loadings: %f',diff_behav_loadings));
            
            % remove the scripts & output directory
            rmdir(out_dir,'s');
            rmpath(fullfile(root_dir,'examples'));
        end
    end
end