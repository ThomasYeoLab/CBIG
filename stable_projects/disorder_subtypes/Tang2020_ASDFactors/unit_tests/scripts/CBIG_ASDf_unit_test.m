 classdef CBIG_ASDf_unit_test < matlab.unittest.TestCase
 % Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
 
     methods (Test)
         function test_short(testCase)
             % create output folder
             CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
             cur_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', ...
                 'disorder_subtypes', 'Tang2020_ASDFactors', 'unit_tests');
             out_dir = [cur_dir '/output'];
             mkdir(out_dir)
             
             % run the example
             cmd = ['sh ' fullfile(cur_dir, 'scripts/', ...
                 'CBIG_ASDf_unit_test.sh') ' ' out_dir];
             system(cmd);
             final_out_file = fullfile(out_dir, 'inference', 'k2r4_factorComp_inf-gamma.dat');
             
             while 1
                 if exist(final_out_file, 'file') == 2
                     break
                 end
             end
             
             % check the results
             ref_dir = fullfile('/mnt', 'eql', 'yeo1', 'CBIG_private_unit_tests_data', ...
                 'stable_projects', 'disorder_subtypes', 'Tang2020_ASDFactors', ...
                 'results', 'results_short');
             curr_out_dir = fullfile(out_dir, 'visualizeFactors', 'k2', 'r4');
             curr_ref_dir = fullfile(ref_dir, 'visualizeFactors', 'k2', 'r4');
             
             % compare mean1 estimate
             mean1 = load(fullfile(curr_out_dir, 'mean1.mat'));
             ref_mean1 = load(fullfile(curr_ref_dir, 'mean1.mat'));
             diff_mean1 = max(max(abs(mean1.mean_corrmat - ref_mean1.mean_corrmat)));
             assert(diff_mean1 < 1e-6,sprintf('maximum mean1 difference: %f',diff_mean1))
             
             % compare mean2 estimate
             mean2 = load(fullfile(curr_out_dir, 'mean2.mat'));
             ref_mean2 = load(fullfile(curr_ref_dir, 'mean2.mat'));
             diff_mean2 = max(max(abs(mean2.mean_corrmat - ref_mean2.mean_corrmat)));
             assert(diff_mean2 < 1e-6,sprintf('maximum mean2 difference: %f',diff_mean2))
             
             % compare factorComp
             factorComp = load(fullfile(curr_out_dir, 'factorComp.txt'));
             ref_factorComp = load(fullfile(curr_ref_dir, 'factorComp.txt'));
             diff_factorComp = max(max(abs(factorComp - ref_factorComp)));
             assert(diff_factorComp < 1e-6,sprintf('maximum factorComp difference: %f',diff_factorComp))
             
             % compare inferred gamma file
             curr_out_dir = fullfile(out_dir, 'inference');
             curr_ref_dir = fullfile(ref_dir, '/inference');
             gamma = load(fullfile(curr_out_dir, 'k2r4_factorComp_inf-gamma.dat'));
             ref_gamma = load(fullfile(curr_ref_dir, 'k2r4_factorComp_inf-gamma.dat'));
             diff_gamma = max(max(abs(gamma - ref_gamma)));
             assert(diff_gamma < 1e-6,sprintf('maximum inferred gamma difference: %f',diff_gamma))
             
             % remove the output directory
             rmdir(out_dir, 's')
         end
     end
 end
