classdef CBIG_VK2019_unit_test < matlab.unittest.TestCase
% Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'disorder_subtypes', 'Kebets2019_TransdiagnosticComponents', 'unit_tests');
            out_dir = fullfile(cur_dir, 'output');
            mkdir(out_dir)
            
            % run the example
            addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'disorder_subtypes', 'Kebets2019_TransdiagnosticComponents', 'examples'));
            
            CBIG_VK2019_example_wrapper(out_dir);
            
            % check the results
            ref_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'disorder_subtypes', 'Kebets2019_TransdiagnosticComponents', 'examples', 'correct_output');
            
            % compare RSFC & behavioral scores & loadings
            ref_Params = load(fullfile(ref_dir, 'PLSresults_example.mat'));          
            ref_Lx = ref_Params.Lx;
            ref_Ly = ref_Params.Ly;
            ref_RSFC_loadings = ref_Params.LC_RSFC_loadings;
            ref_behav_loadings = ref_Params.LC_behav_loadings;
            
            test_Params = load(fullfile(out_dir, 'PLSresults_example.mat'));
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
            rmdir(out_dir, 's');
            rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                'disorder_subtypes', 'Kebets2019_TransdiagnosticComponents', 'examples'));
        end
    end
end