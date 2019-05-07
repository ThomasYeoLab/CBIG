classdef CBIG_MMLDA_unit_test < matlab.unittest.TestCase
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    methods (Test)
        function test_example(testCase)
            % create output folder
            CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
            cur_dir = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/unit_tests'];
            out_dir = [cur_dir '/output'];
            mkdir(out_dir)
            % run the example
            addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/examples'])
            queue = 'circ-spool'
            CBIG_MMLDA_example_wrapper(out_dir, queue)
            % check the results
            ref_dir = [CBIG_CODE_DIR '/stable_projects/disorder_subtypes/' ...
            'Sun2019_ADJointFactors/examples/correct_output'];
            curr_out_dir = [out_dir '/estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3'];
            curr_ref_dir = [ref_dir '/estimation/ADNI2_bl_AD_meanCNstdALL_plus1_2sub/k2/r3'];
            % compare beta1 estimate
            beta1 = load([curr_out_dir '/final.beta1']);
            ref_beta1 = load([curr_ref_dir '/final.beta1']);
            rel_diff_beta1 = max(max(abs((beta1 - ref_beta1) ./ ref_beta1)));
            assert(rel_diff_beta1 < 1e-4,sprintf('maximum beta1 rel difference: %f',rel_diff_beta1))
            % compare beta2 estimate
            beta2 = load([curr_out_dir '/final.beta2']);
            ref_beta2 = load([curr_ref_dir '/final.beta2']);
            rel_diff_beta2 = max(max(abs((beta2 - ref_beta2) ./ ref_beta2)));
            assert(rel_diff_beta2 < 1e-4,sprintf('maximum beta2 rel difference: %f',rel_diff_beta2))
            % compare gamma estimate
            gamma = load([curr_out_dir '/final.gamma']);
            ref_gamma = load([curr_ref_dir '/final.gamma']);
            rel_diff_gamma = max(max(abs((gamma - ref_gamma) ./ ref_gamma)));
            assert(rel_diff_gamma < 1e-4,sprintf('maximum gamma rel difference: %f',rel_diff_gamma))
            % compare likelihood estimate
            likelihood = load([curr_out_dir '/likelihood.dat']);
            likelihood = likelihood(:, 1);
            ref_likelihood = load([curr_ref_dir '/likelihood.dat']);
            ref_likelihood = ref_likelihood(:, 1);
            rel_diff_likelihood = max(max(abs((likelihood - ref_likelihood) ./ ref_likelihood)));
            assert(rel_diff_likelihood < 1e-4,sprintf('maximum likelihood rel difference: %f',rel_diff_likelihood))
            % compare inferred gamma file
            curr_out_dir = [out_dir '/inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];
            curr_ref_dir = [ref_dir '/inference/ADNI2_bl_AD_meanCNstdALL_plus1_2sub'];
            gamma = load([curr_out_dir '/k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat']);
            ref_gamma = load([curr_ref_dir '/k2_inf_ADNI2_bl_AD_meanCNstdALL_plus1_1sub-gamma.dat']);
            rel_diff_gamma = max(max(abs((gamma - ref_gamma) ./ ref_gamma)));
            assert(rel_diff_gamma < 1e-4,sprintf('maximum inferred gamma rel difference: %f',rel_diff_gamma))
            % remove the output directory
            rmdir(out_dir, 's')
            rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/examples'])
        end
    end
end